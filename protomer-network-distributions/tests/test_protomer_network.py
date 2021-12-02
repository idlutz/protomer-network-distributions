import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import tellurium as te
import unittest
from protomer_network_distributions import protomer_network

class Test_Mechanism(unittest.TestCase):
	"""Unit test for protomer_network class"""
	def test___init___(self):
		input_dict = {
						'A_B':10*10**-9,
						'A_C':100*10**-9,
						'B_C':1*10**-9
					 }

	self.test_class = protomer_network(input_dict)
	self.assertTrue(self.test_class.protomer_dict is not None)
	self.assertTrue(self.test_class.simulations is not None)
	self.assertTrue(self.test_class.simulation_distributions is not None)

if __name__ == '__main__':
	unittest.main()
