import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import tellurium as te
import unittest
from protomer_network_distributions import protomer_network

TEST_INPUT_DICT = {
                        'A_B':10*10**-9,
                        'A_C':100*10**-9,
                        'B_C':1*10**-9
                  }

TEST_SUB_DICT = {
                        'A':{},
                        'B':{'label':True},
                        'C':{'label':False}
                }

class Test_VisualizeEquilibriumDistribution(unittest.TestCase):
        """Unit test for VisualizeEquilibriumDistribution method"""
        def test___init___(self):
          self.test_class = protomer_network(TEST_INPUT_DICT)
          self.test_class.SimulateEquilibriumSpeciesDistribution(TEST_SUB_DICT, model_summary=False, plot=False)
          self.test_visualize_output = self.test_class.VisualizeEquilibriumDistribution()
          self.assertTrue(self.test_visualize_output == None)

if __name__ == '__main__':
        unittest.main()
