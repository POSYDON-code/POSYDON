from posydon.binary_evol.singlestar import SingleStar
from posydon.utils.common_functions import check_state_of_star
import unittest

# TODO: this is outdated
Target_state =[
        'H-ZAMS', 'He-ZAMS'
        'Near_H-ZAMS',
        'H-MS',
        'Shell_H_burning',
        'H-rich Core_He_burning', 'stripped_He Core_He_burning',
        'H-rich Shell_He_burning', 'stripped_He Core_He_burning',
        'H-rich Core_C_burning', 'stripped_He Core_C_burning',
        'H-rich Central_C_depletion', 'stripped_He Central_C_depletion'
]

# class Test_check_state_of_star(unittest.TestCase):
#     def test_evolutionary_state(self):
#         attr = {
#             "center_h1": 0.5,
#             "center_he4": 0.48,
#             "center_c12": 0.0,
#             "surface_h1": 0.7,
#             "log_Lnuc": 0.7,
#             "log_LH": 0.7,
#             "c12_c12": 0.0}
#
#         #out = check_state_of_star(star)
#         star = SingleStar(**attr)
#         state = check_state_of_star(star)
#         self.assertIn(state, Target_state)

if __name__ == "__main__":
    unittest.main()
