import unittest
from pyCHSolver import solver
import numpy as np

class TestLaplace(unittest.TestCase):
    def test_random(self):
        a = np.array([[7, 9, 3], [8, 0, 2], [4, 8, 3]]) #random numbers from C code
        a_solution = np.array([[27, -153, 45], [-117, 243, 0], [72, -153, 36]])
        self.assertIsNone(np.testing.assert_array_equal(solver.laplace(a, 3, 3), a_solution))

    def test_identity(self):
        i = np.identity(3) #identity matrix
        i_solution = np.array([[-18, 18, 0], [18, -36, 18], [0, 18, -18]])
        self.assertIsNone(np.testing.assert_array_equal(solver.laplace(i, 3, 3), i_solution))

class TestSource(unittest.TestCase):

    def test_identity_c_old(self):
        i = np.identity(3)  # identity matrix
        sc_solution = np.array([[118, -18, 0], [-18, 136, -18], [0, -18, 118]])
        smu_solution = np.zeros((3,3))
        self.assertIsNone(np.testing.assert_array_equal(solver.source(i, nx = 3, ny = 3, dt = .01), [sc_solution, smu_solution]))

    def test_random_c_old(self):
        a = np.array([[7, 9, 3], [8, 0, 2], [4, 8, 3]]) #random numbers from C code
        sc_solution = np.array([[673, 1053, 255], [917, -243, 200], [328, 953, 264]])
        smu_solution = np.zeros((3,3))

        self.assertIsNone(np.testing.assert_array_equal(solver.source(a, nx = 3, ny = 3, dt = 0.01), [sc_solution, smu_solution]))

    def test_zeros_c_old(self):

        z = np.zeros((3,3)) #zero matrix
        sc_solution = np.zeros((3,3))
        smu_solution = np.zeros((3,3))

        self.assertIsNone(np.testing.assert_array_equal(solver.source(z, nx = 3, ny = 3), [sc_solution, smu_solution]))

class TestRelax(unittest.TestCase):
    def test_identity_c_new(self):
        i = np.identity(3)  # identity matrix
        mu = np.zeros((3,3))
        su = np.array([[118, -18, 0], [-18, 136, -18], [0, -18, 118]])
        sw = np.zeros((3,3))
        c_relax = 2
        xright = 0
        xleft = 1
        dt = 0.01
        Cahn = 0.06**2
        c_new_py, mu_new_py = solver.relax(i, mu, su, sw, 3, 3,
              c_relax=c_relax, xright=xright, xleft=xleft, dt=dt, Cahn=Cahn)

        c_new_c = np.array([[0.983668, 0.022119, -0.011844], [0.022119, 0.965980, 0.020646], [-0.011844, 0.020646, 0.985618]])
        mu_new_c = np.array([[1.019657, -0.060073, 0.000683], [-0.060073, 1.029385, -0.059907], [0.000683, - 0.059907, 1.019993]])
        self.assertIsNone(np.testing.assert_allclose(c_new_py, c_new_c, rtol = 1e-3, atol = 0))
        self.assertIsNone(np.testing.assert_allclose(mu_new_py, mu_new_c, rtol = 1e-3, atol = 0))
class TestnonL(unittest.TestCase):
    def test_identity_c_new(self):
        i = np.identity(4)  # identity matrix
        mu = np.zeros((4,4))
        dt = 0.01
        Cahn = 0.06**2
        ru, rw = solver.nonL(i, mu, 4, 4, dt, Cahn)
        ru_c = 100*np.eye(4)
        rw_c = np.array([[-1.115200,0.115200,0.000000,0.000000],[0.115200,-1.230400,0.115200,0.000000],[0.000000,0.115200,-1.230400,0.115200],[0.000000,0.000000,0.115200,-1.115200]])

        self.assertIsNone(np.testing.assert_array_equal(ru_c, ru))
        self.assertIsNone(np.testing.assert_array_equal(rw, rw_c))

    def test_ones_c_mu(self):
        c = np.array([[1, 1, 1, -1], [-1, -1, -1, -1], [1, 1, -1, 1], [-1, -1, 1, 1]])
        mu = np.array([[1, 1, -1, -1], [1, 1, 1, 1], [-1, 1, 1, -1], [-1, 1, 1, -1]])
        dt = 0.01
        Cahn = 0.06**2
        ru, rw = solver.nonL(c, mu, 4, 4, dt, Cahn)
        ru_c = np.array([[100, 132, 36, -132], [-68, -100, -68, -36], [36, 132, -68, 36], [-132, -68, 132, 68]])

        rw_c = np.array([[-0.115200, -0.115200, -2.230400, 0.115200], [2.230400, 2.230400, 2.115200, 2.115200], [-2.230400, -0.345600, 2.345600, -2.230400], [0.115200, 2.230400, -0.230400, -2.000000]])
        self.assertIsNone(np.testing.assert_array_equal(ru_c, ru))
        self.assertIsNone(np.testing.assert_array_equal(rw, rw_c))

class TestDefect(unittest.TestCase):
    def test_defect(self):
        nx = ny = 4
        nxc = nyc = int(nx/2)
        c = np.array([[1, 1, 1, -1], [-1, -1, -1, -1], [1, 1, -1, 1], [-1, -1, 1, 1]])
        mu = np.array([[1, 1, -1, -1], [1, 1, 1, 1], [-1, 1, 1, -1], [-1, 1, 1, -1]])
        dt = 0.01
        Cahn = 0.06**2
        sc = np.zeros((nx,ny))
        smu = np.zeros((nx,ny))
        uc_new, wc_new = solver.restrict_ch(c, mu, nxc, nyc)
        duc, dwc = solver.defect(c, mu, sc, smu, nx, ny, uc_new, wc_new, nxc, nyc)
        duc_c = np.array([[-8,-4],[4,8]])
        dwc_c = np.array([[-0.064800,-0.382200],[0.064800,0.382200]])
        self.assertIsNone(np.testing.assert_allclose(duc_c, duc, rtol = 1e-10, atol = 0))
        self.assertIsNone(np.testing.assert_allclose(dwc_c, dwc, rtol = 1e-10, atol = 0))

if __name__ == '__main__':
    unittest.main()

#a = np.array([[7, 9, 3, 8, 0, 2, 4, 8], [3, 9, 0, 5, 2, 2, 7, 3], [7, 9, 0, 2, 3, 9, 9, 7], [0, 3, 9, 8, 6, 5, 7, 6], [2, 7, 0, 3, 9, 9, 9, 1], [7, 2, 3, 6, 5, 5, 8, 1], [4, 7, 1, 3, 8, 4, 8, 0], [4, 6, 0, 3, 2, 6, 9, 4]])
