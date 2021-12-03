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

LOWER_BOUND = 6e-09
UPPER_BOUND = 7e-09

class Test_PredictBinaryOutput(unittest.TestCase):
        """Unit test for PredictBinaryOutput method"""
        def test___init___(self):
          self.test_class = protomer_network(TEST_INPUT_DICT)
          self.test_class.SimulateEquilibriumSpeciesDistribution(TEST_SUB_DICT, model_summary=False, plot=False)
          self.test_binary_output1 = self.test_class.PredictBinaryOutput()
          self.assertTrue(self.test_binary_output1 == True)
          self.test_binary_output2 = self.test_class.PredictBinaryOutput(return_cumulative_signal=True)
          self.assertTrue(LOWER_BOUND < self.test_binary_output2 < UPPER_BOUND)

if __name__ == '__main__':
        unittest.main()
