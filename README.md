# First-principles simulation of electron transport and thermoelectric property of materials including electron-phonon scattering, defect scattering, and phonon drag

Jiawei Zhou, Te-Huan Liu, Qichen Song, Qian Xu, Zhiwei Ding, Bolin Liao, and Gang Chen

Department of Mechanical Engineering
Massachusetts Institute of Technology
Cambridge, MA 02139 

# Disclaimer:

This EPW code is a modified version of the EPW v4 code from the open-source Quantum ESPRESSO suite (based on QE 5.4.0), and is released under GNU General Public License (v2). The source codes are published both at Materials Cloud (https://archive.materialscloud.org/record/2020.106, DOI: 10.24435/materialscloud:5a-7s) and GitHub (for comments and questions). The original EPW v4 is developed by S. Poncé, E.R. Margine, C. Verdi, and, F. Giustino, initially released inside Quantum ESPRESSO in 2016. This modified version is dedicated to the simulation of electron-phonon transport properties in quantum materials. Specifically, it calculates the electron-phonon and electron-defect scattering rates and uses them as inputs in Boltzmann transport equation to obtain transport properties (e.g. electrical conductivity, mobility, Seebeck coefficient, thermoelectric power factor, and electronic thermal conductivity).

As part of the Quantum ESPRESSO suite, please cite:

"P. Giannozzi et al., J. Phys.:Condens. Matter 21 395502 (2009); URL http://www.quantum-espresso.org",

in publications or presentations arising from this work. More details at http://www.quantum-espresso.org/quote.

Please also consider citing the EPW papers:

1) F. Giustino, M. L. Cohen, and S. G. Louie, Phys. Rev. B 76, 165108 (2007)
2) S. Poncé, E. R. Margine, C. Verdi, and F. Giustino, Comput. Phys. Comm. 209, 116 (2016)

Our modified EPW code can be cited as:

3) J.W. Zhou, T.-H. Liu, Q.C. Song, Q. Xu, Z.W. Ding, B.L. Liao, and G. Chen, Materials Cloud Archive 2020.106 (2020), doi: 10.24435/materialscloud:5a-7s.

For electron transport property calculations in materials with strong spin-orbit coupling or polar optical phonon scattering, please consider citing:

4) T.-H. Liu, J.W. Zhou, M.D. Li, Z.W. Ding, Q.C. Song, B.L. Liao, L. Fu, and G. Chen, Proceedings of National Academy of Sciences, 115, 879 (2018)
5) T.-H. Liu, J.W. Zhou, B.L. Liao, D.J. Singh, and G. Chen, Phys. Rev. B, 95, 075206 (2017)

For phonon drag calculations, please consider citing:

6) J.W. Zhou, B.L. Liao, B. Qiu, S. Huberman, K. Esfarjani, M.S. Dresselhaus and G. Chen, Proceedings of National Academy of Sciences, 112, 14777 (2015)

For questions related to the general usage of EPW and the functions of the original EPW code, please refer to the EPW website (https://epw-code.org/) and references there. For questions related to the modified EPW code, please address to Issues at our GitHub repository page. Comments to the code are also welcome.

Acknowledgement: The codes were developed over the years under the support of the following programs: S3TEC, an Energy Frontier Research Center funded by the U.S. Department of Energy, Office of Basic Energy Sciences under Award No. DE- SC0001299/DE-FG02-09ER46577 and the DARPA MATRIX program, under Grant HR0011-16-2-0041.

