#ifndef GLOBAL_FP_H
#define GLOBAL_FP_H

#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_lcao/module_gint/grid_technique.h"
#include "local_orbital_wfc.h"
#include "local_orbital_charge.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_gen_fixedH.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_hamilt.h" 
#include "module_basis/module_ao/ORB_read.h"
#include "module_basis/module_ao/ORB_gen_tables.h"
#ifdef __EXX
#include "module_ri/Exx_LRI.h"
#endif

namespace GlobalC
{
extern Grid_Driver GridD;
#ifdef __EXX
extern Exx_LRI<double> exx_lri_double; // Peize Lin add 2022-08-06
extern Exx_LRI<std::complex<double>> exx_lri_complex; // Peize Lin add 2022-08-06
#endif
}
#endif
