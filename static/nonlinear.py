import fenics as fe
import sys
import argparse
import os
import datetime
import numpy as np
import sys
from pathlib import Path
from tqdm import tqdm
import matplotlib.pyplot as plt
import toml

fe.set_log_level(fe.LogLevel.ERROR)
fe.parameters["reorder_dofs_serial"] = False

parser = argparse.ArgumentParser()

parser.add_argument("--opt", choices=["output", "debug"], default="output")
parser.add_argument("--machine", type=str, default= 'data/machine/general.toml')
parser.add_argument("--material", type=str, default= 'data/material/beam.toml')
parser.add_argument("--model", type=str, default= 'data/model/beamShort.toml')
parser.add_argument("--boundary", type=str, default= 'data/boundary/topology_static.toml')
parser.add_argument("--mode", type=str, choices=["optimization", "check"], default="optimization")

args = parser.parse_args()


class topology() :
    def __init__(self, machine, material, model, boundary) :
        self.dir_mesh = machine["dir_mesh"]
        self.dir_consequence = machine["dir_consequence"]
        if EXPORT :
            date = datetime.datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
            program_name = os.path.splitext(os.path.basename(__file__))[0]
            folder_name = Path(program_name) / Path(date)
            self.dir_consequence = self.dir_consequence / folder_name
            os.makedirs(self.dir_consequence, exist_ok=True)
        self.mesh_name = model["mesh_name"]
        self.region_name = model["region_name"]
        self.label_fix = model["label_fix"]
        self.label_ext = model["label_ext"]
        self.dim = model["dim"]

        self.set_material_info(material)
        self.set_boundary_info(boundary)
        self.output_times = 0
        self.obj_0_history = np.zeros(dtype=np.float64, shape=0)
        self.obj_vol = 0
        self.psi_max_history = np.zeros(dtype=np.float64, shape=0)
        
        if EXPORT :
            self.file_results = fe.XDMFFile(str(self.dir_consequence / Path("SOLID.xdmf")))
            self.file_results.parameters["flush_output"] = True
            self.file_results.parameters["functions_share_mesh"] = True
    
    def set_material_info(self, material) :
        self.young = material["young"]
        self.nu = material["nu"]
        self.mu = self.young / (2 * (1 - 2 * self.nu))
        self.la = self.young * self.nu / ((1 + self.nu) * (1 - 2 * self.nu))

    def set_boundary_info(self, boundary) :
        structure = boundary["structure"]
        optimization = boundary["optimization"]
        self.traction_force = structure["traction_force"]
        self.alpha_exp = optimization["alpha_exp"]
        self.alpha_vol = optimization["alpha_vol"]
        self.ca = optimization["ca"]
        self.cD = optimization["cD"]
        self.max_number_iter = optimization["max_number_iter"]
        self.max_move = optimization["max_move"]
        self.output_span = optimization["output_span"]
        self.plot_output_span = optimization["plot_output_span"]

    def define_mesh(self) :
        self.mesh = fe.Mesh(str(self.dir_mesh / Path(self.mesh_name)))
        self.facet_function = fe.MeshFunction("size_t", self.mesh, str(self.dir_mesh / Path(self.region_name)))
        self.ds = fe.Measure("ds", domain=self.mesh, subdomain_data=self.facet_function)
        self.dx = fe.dx(metadata={"quadrature_degree": 6})
        
    def set_mesh_function(self) :
        self.vector_space = fe.VectorFunctionSpace(self.mesh, "CG", 2)
        self.scaler_space = fe.FunctionSpace(self.mesh, "CG", 1)
        self.dis = fe.Function(self.vector_space, name="displacement")
        self.delta_dis = fe.Function(self.vector_space, name="delta_dis")
        self.vec_test = fe.TestFunction(self.vector_space)
        self.vec_trial = fe.TrialFunction(self.vector_space)
        self.theta = fe.Function(self.scaler_space, name="theta")
        self.phi = fe.Function(self.scaler_space, name="phi")  
        self.psi = fe.Function(self.scaler_space, name="psi")  
        self.psi_0 = fe.Function(self.scaler_space, name="psi_0")  
        self.psi_vol = fe.Function(self.scaler_space, name="psi_vol")  
        self.scl_test = fe.TestFunction(self.scaler_space)
        self.scl_trial = fe.TrialFunction(self.scaler_space)
        self.theta.interpolate(fe.Constant(0.0))
        
    def set_boundary(self) :
        self.traction_force = fe.Constant((0.0, self.traction_force))
        self.boundary_dis = fe.Constant((0.0, 0.0))
        self.bc = fe.DirichletBC(self.vector_space, self.boundary_dis, self.facet_function, self.label_fix)

    def get_lagrangean(self, theta, dis, delta_dis) :
        obj0 = self.get_obj0(theta, dis)
        R = self.get_R(theta, dis, delta_dis)
        L = obj0 - R
        return L

    def get_obj0(self, theta, dis) :
        return fe.inner(self.traction_force, dis) * self.ds(self.label_ext)

    def get_rho_hat(self, theta) :
        return self.get_phi(theta)**self.alpha_exp
    
    def get_F(self, u) :
        return fe.Identity(self.dim) + fe.grad(u)
    
    def get_P_neohook(self, u) :
        F = self.get_F(u)
        J = fe.det(F)
        P = self.mu * F + (- self.mu + self.la * fe.ln(J)) * fe.inv(F).T
        return P
    
    def get_dF(self, du) :
        return fe.grad(du)

    def get_inner_term(self, theta, dis, delta_dis) :
        return self.get_rho_hat(theta) * fe.inner(self.get_P_neohook(dis), self.get_dF(delta_dis)) * self.dx
    
    def get_exter_term(self, theta, dis, delta_dis) :
        return fe.inner(self.traction_force, delta_dis) * self.ds(self.label_ext)
        
    def get_R(self, theta, dis, delta_dis) :
        inner_term = self.get_inner_term(theta, dis, delta_dis)
        exter_term = self.get_exter_term(theta, dis, delta_dis)
        return - inner_term + exter_term

    def solve_structure_problem(self) :
        # delta = self.get_R(self.theta, self.dis, self.vec_test)
        L = self.get_lagrangean(self.theta, self.dis, self.delta_dis)
        delta = fe.derivative(L, self.delta_dis, self.vec_test)
        Hesse = fe.derivative(delta, self.dis, self.vec_trial)
        problem = fe.NonlinearVariationalProblem(delta, self.dis, self.bc, J= Hesse)
        solver = fe.NonlinearVariationalSolver(problem)
        prm = solver.parameters
        solver.solve()

    def solve_adjoint_problem(self) :
        L = self.get_lagrangean(self.theta, self.dis, self.delta_dis)
        delta = fe.derivative(L, self.dis, self.vec_test)
        Hesse = fe.derivative(delta, self.delta_dis, self.vec_trial)
        problem = fe.NonlinearVariationalProblem(delta, self.delta_dis, self.bc, J= Hesse)
        solver = fe.NonlinearVariationalSolver(problem)
        prm = solver.parameters
        solver.solve()

    def get_sens_0(self) :
        L = self.get_lagrangean(self.theta, self.dis, self.delta_dis)
        sens0 = fe.derivative(L, self.theta, self.scl_test)
        return sens0

    def get_sens_vol(self, theta, vartheta) :
        L = self.get_phi(theta) * self.dx
        sens_vol = fe.derivative(L, theta, vartheta)
        return sens_vol

    def get_a_h1(self, scl_trial, scl_test) :
        return self.ca * (fe.inner(fe.grad(scl_trial), fe.grad(scl_test)) + self.cD * scl_trial * scl_test) * self.dx
        
    def cal_psi_0(self) :
        a = self.get_a_h1(self.scl_trial, self.scl_test)
        l = - self.get_sens_0()
        fe.solve(a == l, self.psi_0)
        
    def cal_psi_vol(self) :
        a = self.get_a_h1(self.scl_trial, self.scl_test)
        l = - self.get_sens_vol(self.theta, self.scl_test)
        fe.solve(a == l, self.psi_vol)
        
    def cal_phi(self):
        a = self.scl_trial * self.scl_test * fe.dx
        l = self.get_phi(self.theta)**self.alpha_exp * self.scl_test * fe.dx
        fe.solve(a == l, self.phi)

    def get_product(self, theta, vartheta) :
        return fe.assemble(self.get_sens_vol(theta, vartheta))
        
    def update_theta(self) :
        product0 = self.get_product(self.theta, self.psi_0)
        product1 = self.get_product(self.theta, self.psi_vol)
        lambda1 = - (self.obj_vol + product0) / product1
        self.psi.vector()[:] = self.psi_0.vector()[:] + lambda1 * self.psi_vol.vector()[:]
        psi_max = np.max(np.abs(self.psi.vector().get_local()))
        self.theta.vector().axpy(1.0, self.psi.vector())
        
        if EXPORT :
            self.psi_max_history = np.append(self.psi_max_history, psi_max)
            plt.cla()
            plt.plot(self.psi_max_history)
            plt.savefig(self.dir_consequence / Path("psi_max"))
        
    def get_sinh(self, x) :
        return (fe.exp(x) - fe.exp(- x)) / 2

    def get_cosh(self, x) :
        return (fe.exp(x) + fe.exp(- x)) / 2
        
    def get_tanh(self, x) :
        return self.get_sinh(x) / self.get_cosh(x)
        
    def get_phi(self, theta) :
        return 0.5 * self.get_tanh(theta) + 0.5

    def cal_vol_crn(self) :
        return fe.assemble(self.get_phi(self.theta) * fe.dx)
    
    def get_obj_vol(self) :
        return self.cal_vol_crn() - self.vol_init

    def cal_obj_vol(self) :
        self.obj_vol = self.cal_vol_crn() - self.vol_init
    
    def push_xdmf_data(self) :
        if EXPORT and self.time_step % self.output_span == 0 :
            self.file_results.write(self.dis, self.output_times)
            self.file_results.write(self.theta, self.output_times)
            self.file_results.write(self.phi, self.output_times)
            self.output_times += 1
    
    def plot_phi_field(self) :
        if EXPORT and self.time_step % self.plot_output_span == 0 :
            plt.figure(figsize=(10, 8))
            plt.clf()
            
            values = self.phi.compute_vertex_values(self.mesh)
            coordinates = self.mesh.coordinates()
            cells = self.mesh.cells()
            
            c = plt.tripcolor(
                coordinates[:, 0],
                coordinates[:, 1],
                cells,
                values,
                cmap="coolwarm",
                shading="gouraud",
            )
            
            plt.colorbar(c, label="Phi")
            plt.title(f"Phi Field (Iteration {self.time_step})")
            plt.axis("equal")
            plt.xlabel("X")
            plt.ylabel("Y")
            plt.savefig(self.dir_consequence / Path("phi_field"))
            plt.close()

    def cal_obj_0(self) :
        self.obj0 = self.get_obj0(self.theta, self.dis) 
        self.obj0 = fe.assemble(self.obj0)
        self.obj_0_history = np.append(self.obj_0_history, self.obj0)
        if EXPORT :
            plt.cla()
            plt.plot(self.obj_0_history)
            plt.savefig(self.dir_consequence / Path("object_history"))
        
    def main_optimization(self) :
        self.define_mesh()
        self.set_mesh_function()
        self.set_boundary()
        
        self.vol_init = self.cal_vol_crn()
        
        for _ in tqdm(range(self.max_number_iter)) :
            self.time_step = _
            self.solve_structure_problem()
            self.solve_adjoint_problem()
            self.cal_obj_0()
            self.cal_obj_vol()
            self.cal_psi_0()
            self.cal_psi_vol()
            self.update_theta()
            self.cal_phi()
            self.push_xdmf_data()
            self.plot_phi_field()
            
    def set_init_check_sens(self) :
        self.delta_theta = 1.0e-5
        self.theta_np = self.theta.vector().get_local()
        self.num_p_scl = self.theta_np.shape[0]
        self.num_sample = 10
        self.p_random = np.random.choice(self.num_p_scl, size=self.num_sample, replace=False)
        
    def get_obj_posi_or_nega(self, p, sign) :
        self.theta_np.fill(0)
        self.theta_np[p] += sign * self.delta_theta
        self.theta.vector().set_local(self.theta_np)
        self.solve_structure_problem()
        obj0 = self.get_obj0(self.theta, self.dis)
        obj0 = fe.assemble(obj0)
        return obj0
            
    def cal_sens_from_diff(self) :
        self.sens_from_diff = np.zeros(dtype=np.float64, shape=self.num_p_scl)
        for _p in tqdm(range(self.num_sample)) :
            p = self.p_random[_p]
            obj_posi = self.get_obj_posi_or_nega(p, 1)
            obj_nega = self.get_obj_posi_or_nega(p, -1)
            dfds = (obj_posi - obj_nega) / (2 * self.delta_theta)
            self.sens_from_diff[p] = dfds
            
    def cal_sens_from_deri(self) :
        self.theta.assign(fe.Constant(0.0))
        self.solve_structure_problem()
        self.solve_adjoint_problem()
        self.sens_from_deri = self.get_sens_0()
        self.sens_from_deri = fe.assemble(self.sens_from_deri)
        self.sens_0_np = self.sens_from_deri.get_local()
            
    def main_check_sens(self) :
        self.define_mesh()
        self.set_mesh_function()
        self.set_boundary()
        
        self.vol_init = self.cal_vol_crn()
        
        self.set_init_check_sens()
        self.cal_sens_from_diff()
        self.cal_sens_from_deri()
        
        plt.plot(self.sens_from_diff[self.p_random], label="diff", marker='o', markerfacecolor='none')
        plt.plot(self.sens_0_np[self.p_random], label="sens")
        plt.legend()
        plt.show()
        plt.savefig("sens_diff_nonlinear")
        sys.exit()
        
if __name__ == "__main__" :
    dir_material = Path(args.material)
    dir_machine = Path(args.machine)
    dir_model = Path(args.model)
    dir_boundary = Path(args.boundary)
    mode = args.mode
    EXPORT = args.opt == "output"

    if EXPORT :
        DATE = datetime.datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
        FILE_NAME = os.path.splitext(os.path.basename(__file__))[0]
        # self.dir_consequence =  + FILE_NAME + "/" + DATE
        # self.dir_consequence = "C:/Users/x1x4x/"
        # os.makedirs(self.dir_consequence, exist_ok=True)

    with open(dir_machine, "r") as file :
        machine = toml.load(file)

    with open(dir_material, "r") as file :
        material = toml.load(file)

    with open(dir_model, "r") as file :
        model = toml.load(file)
    
    with open(dir_boundary, "r") as file :
        boundary = toml.load(file)

    obj = topology(machine, material, model, boundary)

    if mode == "optimization" :
        obj.main_optimization()
    elif mode == "check" :
        obj.main_check_sens()
    else :
        sys.exit("error : mode is optimization or check")
    
        

