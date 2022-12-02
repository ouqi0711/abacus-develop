#include "forces.h"

#include "../module_symmetry/symmetry.h"
#include "global.h"
// new
#include "module_base/math_integral.h"
#include "module_base/timer.h"
#include "module_base/tool_threading.h"
#include "module_elecstate/potentials/efield.h"
#include "module_surchem/surchem.h"
#include "module_elecstate/potentials/gatefield.h"
#include "module_vdw/vdw.h"

#ifdef _OPENMP
#include <omp.h>
#endif

double Forces::output_acc = 1.0e-8; // (Ryd/angstrom).

Forces::Forces()
{
}

Forces::~Forces()
{
}

#include "../module_base/mathzone.h"
void Forces::init(ModuleBase::matrix& force, const ModuleBase::matrix& wg, const Charge* const chr, const psi::Psi<std::complex<double>>* psi_in)
{
    ModuleBase::TITLE("Forces", "init");
    this->nat = GlobalC::ucell.nat;
    force.create(nat, 3);

    ModuleBase::matrix forcelc(nat, 3);
    ModuleBase::matrix forceion(nat, 3);
    ModuleBase::matrix forcecc(nat, 3);
    ModuleBase::matrix forcenl(nat, 3);
    ModuleBase::matrix forcescc(nat, 3);
    this->cal_force_loc(forcelc, GlobalC::rhopw, chr);
    this->cal_force_ew(forceion, GlobalC::rhopw);
    this->cal_force_nl(forcenl, wg, psi_in);
    this->cal_force_cc(forcecc, GlobalC::rhopw, chr);
    this->cal_force_scc(forcescc, GlobalC::rhopw);

    ModuleBase::matrix stress_vdw_pw; //.create(3,3);
    ModuleBase::matrix force_vdw;
    force_vdw.create(nat, 3);
    auto vdw_solver = vdw::make_vdw(GlobalC::ucell, INPUT);
    if (vdw_solver != nullptr)
    {
        const std::vector<ModuleBase::Vector3<double>> &force_vdw_temp = vdw_solver->get_force();
        for (int iat = 0; iat < GlobalC::ucell.nat; ++iat)
        {
            force_vdw(iat, 0) = force_vdw_temp[iat].x;
            force_vdw(iat, 1) = force_vdw_temp[iat].y;
            force_vdw(iat, 2) = force_vdw_temp[iat].z;
        }
        if (GlobalV::TEST_FORCE)
        {
            Forces::print("VDW      FORCE (Ry/Bohr)", force_vdw);
        }
    }

    ModuleBase::matrix force_e;
    if (GlobalV::EFIELD_FLAG)
    {
        force_e.create(GlobalC::ucell.nat, 3);
        elecstate::Efield::compute_force(GlobalC::ucell, force_e);
        if (GlobalV::TEST_FORCE)
        {
            Forces::print("EFIELD      FORCE (Ry/Bohr)", force_e);
        }
    }

    ModuleBase::matrix force_gate;
    if (GlobalV::GATE_FLAG)
    {
        force_gate.create(GlobalC::ucell.nat, 3);
        elecstate::Gatefield::compute_force(GlobalC::ucell, force_gate);
        if (GlobalV::TEST_FORCE)
        {
            Forces::print("GATEFIELD      FORCE (Ry/Bohr)", force_gate);
        }
    }

    ModuleBase::matrix forcesol;
    if (GlobalV::imp_sol)
    {
        forcesol.create(GlobalC::ucell.nat, 3);
        GlobalC::solvent_model.cal_force_sol(GlobalC::ucell, GlobalC::rhopw, forcesol);
        if(GlobalV::TEST_FORCE)
        {
            Forces::print("IMP_SOL      FORCE (Ry/Bohr)", forcesol);
        }
    }

    // impose total force = 0
    int iat = 0;
    for (int ipol = 0; ipol < 3; ipol++)
    {
        double sum = 0.0;
        iat = 0;

        for (int it = 0; it < GlobalC::ucell.ntype; it++)
        {
            for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
            {
                force(iat, ipol) = forcelc(iat, ipol) + forceion(iat, ipol) + forcenl(iat, ipol) + forcecc(iat, ipol)
                                   + forcescc(iat, ipol);

                if (vdw_solver != nullptr) // linpz and jiyy added vdw force, modified by zhengdy
                {
                    force(iat, ipol) += force_vdw(iat, ipol);
                }

                if (GlobalV::EFIELD_FLAG)
                {
                    force(iat, ipol) = force(iat, ipol) + force_e(iat, ipol);
                }

                if (GlobalV::GATE_FLAG)
                {
                    force(iat, ipol) = force(iat, ipol) + force_gate(iat, ipol);
                }

                if(GlobalV::imp_sol)
                {
                    force(iat,ipol) = force(iat, ipol) + forcesol(iat, ipol);
                }

				sum += force(iat, ipol);

                iat++;
            }
        }

        if(!(GlobalV::GATE_FLAG || GlobalV::EFIELD_FLAG))
        {
            double compen = sum / GlobalC::ucell.nat;
            for (int iat = 0; iat < GlobalC::ucell.nat; ++iat)
            {
                force(iat, ipol) = force(iat, ipol) - compen;
            }
        }
    }

    if(GlobalV::GATE_FLAG || GlobalV::EFIELD_FLAG)
    {
        GlobalV::ofs_running << "Atomic forces are not shifted if gate_flag or efield_flag == true!" << std::endl;
    }

    if (ModuleSymmetry::Symmetry::symm_flag)
    {
        double* pos;
        double d1, d2, d3;
        pos = new double[GlobalC::ucell.nat * 3];
        ModuleBase::GlobalFunc::ZEROS(pos, GlobalC::ucell.nat * 3);
        int iat = 0;
        for (int it = 0; it < GlobalC::ucell.ntype; it++)
        {
            // Atom* atom = &GlobalC::ucell.atoms[it];
            for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
            {
                pos[3 * iat] = GlobalC::ucell.atoms[it].taud[ia].x;
                pos[3 * iat + 1] = GlobalC::ucell.atoms[it].taud[ia].y;
                pos[3 * iat + 2] = GlobalC::ucell.atoms[it].taud[ia].z;
                for (int k = 0; k < 3; ++k)
                {
                    GlobalC::symm.check_translation(pos[iat * 3 + k], -floor(pos[iat * 3 + k]));
                    GlobalC::symm.check_boundary(pos[iat * 3 + k]);
                }
                iat++;
            }
        }

        for (int iat = 0; iat < GlobalC::ucell.nat; iat++)
        {
            ModuleBase::Mathzone::Cartesian_to_Direct(force(iat, 0),
                                                      force(iat, 1),
                                                      force(iat, 2),
                                                      GlobalC::ucell.a1.x,
                                                      GlobalC::ucell.a1.y,
                                                      GlobalC::ucell.a1.z,
                                                      GlobalC::ucell.a2.x,
                                                      GlobalC::ucell.a2.y,
                                                      GlobalC::ucell.a2.z,
                                                      GlobalC::ucell.a3.x,
                                                      GlobalC::ucell.a3.y,
                                                      GlobalC::ucell.a3.z,
                                                      d1,
                                                      d2,
                                                      d3);

            force(iat, 0) = d1;
            force(iat, 1) = d2;
            force(iat, 2) = d3;
        }
        GlobalC::symm.force_symmetry(force, pos, GlobalC::ucell);
        for (int iat = 0; iat < GlobalC::ucell.nat; iat++)
        {
            ModuleBase::Mathzone::Direct_to_Cartesian(force(iat, 0),
                                                      force(iat, 1),
                                                      force(iat, 2),
                                                      GlobalC::ucell.a1.x,
                                                      GlobalC::ucell.a1.y,
                                                      GlobalC::ucell.a1.z,
                                                      GlobalC::ucell.a2.x,
                                                      GlobalC::ucell.a2.y,
                                                      GlobalC::ucell.a2.z,
                                                      GlobalC::ucell.a3.x,
                                                      GlobalC::ucell.a3.y,
                                                      GlobalC::ucell.a3.z,
                                                      d1,
                                                      d2,
                                                      d3);
            force(iat, 0) = d1;
            force(iat, 1) = d2;
            force(iat, 2) = d3;
        }
        // std::cout << "nrotk =" << GlobalC::symm.nrotk << std::endl;
        delete[] pos;
    }

    GlobalV::ofs_running << std::setiosflags(ios::fixed) << std::setprecision(6) << std::endl;
    /*if(GlobalV::TEST_FORCE)
    {
        Forces::print("LOCAL    FORCE (Ry/Bohr)", forcelc);
        Forces::print("NONLOCAL FORCE (Ry/Bohr)", forcenl);
        Forces::print("NLCC     FORCE (Ry/Bohr)", forcecc);
        Forces::print("ION      FORCE (Ry/Bohr)", forceion);
        Forces::print("SCC      FORCE (Ry/Bohr)", forcescc);
        if(GlobalV::EFIELD) Forces::print("EFIELD   FORCE (Ry/Bohr)", force_e);
    }*/

    /*
        Forces::print("   TOTAL-FORCE (Ry/Bohr)", force);

        if(INPUT.out_force)                                                   // pengfei 2016-12-20
        {
            std::ofstream ofs("FORCE.dat");
            if(!ofs)
            {
                std::cout << "open FORCE.dat error !" <<std::endl;
            }
            for(int iat=0; iat<GlobalC::ucell.nat; iat++)
            {
                ofs << "   " << force(iat,0)*ModuleBase::Ry_to_eV / 0.529177
                    << "   " << force(iat,1)*ModuleBase::Ry_to_eV / 0.529177
                    << "   " << force(iat,2)*ModuleBase::Ry_to_eV / 0.529177 << std::endl;
            }
            ofs.close();
        }
    */

    // output force in unit eV/Angstrom
    GlobalV::ofs_running << std::endl;

	if(GlobalV::TEST_FORCE)
	{
		Forces::print("LOCAL    FORCE (eV/Angstrom)", forcelc,0);
		Forces::print("NONLOCAL FORCE (eV/Angstrom)", forcenl,0);
		Forces::print("NLCC     FORCE (eV/Angstrom)", forcecc,0);
		Forces::print("ION      FORCE (eV/Angstrom)", forceion,0);
		Forces::print("SCC      FORCE (eV/Angstrom)", forcescc,0);
		if(GlobalV::EFIELD_FLAG) Forces::print("EFIELD   FORCE (eV/Angstrom)", force_e,0);
        if(GlobalV::GATE_FLAG) Forces::print("GATEFIELD   FORCE (eV/Angstrom)", force_gate,0);
        if(GlobalV::imp_sol) Forces::print("IMP_SOL   FORCE (eV/Angstrom)", forcesol,0);
	}
	Forces::print("   TOTAL-FORCE (eV/Angstrom)", force,0);

    return;
}

void Forces::print_to_files(std::ofstream& ofs, const std::string& name, const ModuleBase::matrix& f)
{
    int iat = 0;
    ofs << " " << name;
    ofs << std::setprecision(8);
    // ofs << std::setiosflags(ios::showpos);

    double fac = ModuleBase::Ry_to_eV / 0.529177; // (eV/A)

    if (GlobalV::TEST_FORCE)
    {
        std::cout << std::setiosflags(ios::showpos);
        std::cout << " " << name;
        std::cout << std::setprecision(8);
    }

    for (int it = 0; it < GlobalC::ucell.ntype; it++)
    {
        for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
        {
            ofs << " " << std::setw(5) << it << std::setw(8) << ia + 1 << std::setw(20) << f(iat, 0) * fac
                << std::setw(20) << f(iat, 1) * fac << std::setw(20) << f(iat, 2) * fac << std::endl;

            if (GlobalV::TEST_FORCE)
            {
                std::cout << " " << std::setw(5) << it << std::setw(8) << ia + 1 << std::setw(20) << f(iat, 0) * fac
                          << std::setw(20) << f(iat, 1) * fac << std::setw(20) << f(iat, 2) * fac << std::endl;
            }
            iat++;
        }
    }

    GlobalV::ofs_running << std::resetiosflags(ios::showpos);
    std::cout << std::resetiosflags(ios::showpos);
    return;
}

void Forces::print(const std::string& name, const ModuleBase::matrix& f, bool ry)
{
    ModuleBase::GlobalFunc::NEW_PART(name);

    GlobalV::ofs_running << " " << std::setw(8) << "atom" << std::setw(15) << "x" << std::setw(15) << "y"
                         << std::setw(15) << "z" << std::endl;
    GlobalV::ofs_running << std::setiosflags(ios::showpos);
    GlobalV::ofs_running << std::setprecision(8);

    const double fac = ModuleBase::Ry_to_eV / 0.529177;

    if (GlobalV::TEST_FORCE)
    {
        std::cout << " --------------- " << name << " ---------------" << std::endl;
        std::cout << " " << std::setw(8) << "atom" << std::setw(15) << "x" << std::setw(15) << "y" << std::setw(15)
                  << "z" << std::endl;
        std::cout << std::setiosflags(ios::showpos);
        std::cout << std::setprecision(6);
    }

    int iat = 0;
    for (int it = 0; it < GlobalC::ucell.ntype; it++)
    {
        for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
        {
            std::stringstream ss;
            ss << GlobalC::ucell.atoms[it].label << ia + 1;

            if (ry) // output Rydberg Unit
            {
                GlobalV::ofs_running << " " << std::setw(8) << ss.str();
                if (abs(f(iat, 0)) > Forces::output_acc)
                    GlobalV::ofs_running << std::setw(15) << f(iat, 0);
                else
                    GlobalV::ofs_running << std::setw(15) << "0";
                if (abs(f(iat, 1)) > Forces::output_acc)
                    GlobalV::ofs_running << std::setw(15) << f(iat, 1);
                else
                    GlobalV::ofs_running << std::setw(15) << "0";
                if (abs(f(iat, 2)) > Forces::output_acc)
                    GlobalV::ofs_running << std::setw(15) << f(iat, 2);
                else
                    GlobalV::ofs_running << std::setw(15) << "0";
                GlobalV::ofs_running << std::endl;
            }
            else
            {
                GlobalV::ofs_running << " " << std::setw(8) << ss.str();
                if (abs(f(iat, 0)) > Forces::output_acc)
                    GlobalV::ofs_running << std::setw(15) << f(iat, 0) * fac;
                else
                    GlobalV::ofs_running << std::setw(15) << "0";
                if (abs(f(iat, 1)) > Forces::output_acc)
                    GlobalV::ofs_running << std::setw(15) << f(iat, 1) * fac;
                else
                    GlobalV::ofs_running << std::setw(15) << "0";
                if (abs(f(iat, 2)) > Forces::output_acc)
                    GlobalV::ofs_running << std::setw(15) << f(iat, 2) * fac;
                else
                    GlobalV::ofs_running << std::setw(15) << "0";
                GlobalV::ofs_running << std::endl;
            }

            if (GlobalV::TEST_FORCE && ry)
            {
                std::cout << " " << std::setw(8) << ss.str();
                if (abs(f(iat, 0)) > Forces::output_acc)
                    std::cout << std::setw(15) << f(iat, 0);
                else
                    std::cout << std::setw(15) << "0";
                if (abs(f(iat, 1)) > Forces::output_acc)
                    std::cout << std::setw(15) << f(iat, 1);
                else
                    std::cout << std::setw(15) << "0";
                if (abs(f(iat, 2)) > Forces::output_acc)
                    std::cout << std::setw(15) << f(iat, 2);
                else
                    std::cout << std::setw(15) << "0";
                std::cout << std::endl;
            }
            else if (GlobalV::TEST_FORCE)
            {
                std::cout << " " << std::setw(8) << ss.str();
                if (abs(f(iat, 0)) > Forces::output_acc)
                    std::cout << std::setw(15) << f(iat, 0) * fac;
                else
                    std::cout << std::setw(15) << "0";
                if (abs(f(iat, 1)) > Forces::output_acc)
                    std::cout << std::setw(15) << f(iat, 1) * fac;
                else
                    std::cout << std::setw(15) << "0";
                if (abs(f(iat, 2)) > Forces::output_acc)
                    std::cout << std::setw(15) << f(iat, 2) * fac;
                else
                    std::cout << std::setw(15) << "0";
                std::cout << std::endl;
            }

            iat++;
        }
    }

    GlobalV::ofs_running << std::resetiosflags(ios::showpos);
    std::cout << std::resetiosflags(ios::showpos);
    return;
}

void Forces::cal_force_loc(ModuleBase::matrix& forcelc, ModulePW::PW_Basis* rho_basis, const Charge* const chr)
{
    ModuleBase::TITLE("Forces", "cal_force_loc");
    ModuleBase::timer::tick("Forces", "cal_force_loc");

    std::complex<double>* aux = new std::complex<double>[rho_basis->nmaxgr];
    // now, in all pools , the charge are the same,
    // so, the force calculated by each pool is equal.


    /*
        blocking rho_basis->nrxx for data locality.

        By blocking aux with block size 1024,
        we can keep the blocked aux in L1 cache when iterating GlobalV::NSPIN loop
        performance will be better when number of atom is quite huge
    */
    const int block_ir = 1024;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int irb = 0; irb < rho_basis->nrxx; irb += block_ir)
    {
        // calculate the actual task length of this block
        int ir_end = std::min(irb + block_ir, rho_basis->nrxx);

        { // is = 0
            for (int ir = irb; ir < ir_end; ++ir)
            { // initialize aux
                aux[ir] = std::complex<double>(chr->rho[0][ir], 0.0);
            }
        }
        for (int is = 1; is < GlobalV::NSPIN; is++)
        {
            for (int ir = irb; ir < ir_end; ++ir)
            { // accumulate aux
                aux[ir] += std::complex<double>(chr->rho[is][ir], 0.0);
            }
        }
    }

    // to G space. maybe need fftw with OpenMP
    rho_basis->real2recip(aux, aux);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int iat = 0; iat < GlobalC::ucell.nat; ++iat)
    {
        // read `it` `ia` from the table
        int it = GlobalC::ucell.iat2it[iat];
        int ia = GlobalC::ucell.iat2ia[iat];
        for (int ig = 0; ig < rho_basis->npw; ig++)
        {
            const double phase = ModuleBase::TWO_PI * (rho_basis->gcar[ig] * GlobalC::ucell.atoms[it].tau[ia]);
            const double factor = GlobalC::ppcell.vloc(it, rho_basis->ig2igg[ig])
                                  * (cos(phase) * aux[ig].imag() + sin(phase) * aux[ig].real());
            forcelc(iat, 0) += rho_basis->gcar[ig][0] * factor;
            forcelc(iat, 1) += rho_basis->gcar[ig][1] * factor;
            forcelc(iat, 2) += rho_basis->gcar[ig][2] * factor;
        }
        forcelc(iat, 0) *= (GlobalC::ucell.tpiba * GlobalC::ucell.omega);
        forcelc(iat, 1) *= (GlobalC::ucell.tpiba * GlobalC::ucell.omega);
        forcelc(iat, 2) *= (GlobalC::ucell.tpiba * GlobalC::ucell.omega);
    }

    // this->print(GlobalV::ofs_running, "local forces", forcelc);
    Parallel_Reduce::reduce_double_pool(forcelc.c, forcelc.nr * forcelc.nc);
    delete[] aux;
    ModuleBase::timer::tick("Forces", "cal_force_loc");
    return;
}

#include "H_Ewald_pw.h"
void Forces::cal_force_ew(ModuleBase::matrix& forceion, ModulePW::PW_Basis* rho_basis)
{
    ModuleBase::TITLE("Forces", "cal_force_ew");
    ModuleBase::timer::tick("Forces", "cal_force_ew");

    double fact = 2.0;
    std::complex<double>* aux = new std::complex<double>[rho_basis->npw];

    /*
        blocking rho_basis->nrxnpwx for data locality.

        By blocking aux with block size 1024,
        we can keep the blocked aux in L1 cache when iterating GlobalC::ucell.ntype loop
        performance will be better when number of atom is quite huge
    */
    const int block_ig = 1024;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int igb = 0; igb < rho_basis->npw; igb += block_ig)
    {
        // calculate the actual task length of this block
        int ig_end = std::min(igb + block_ig, rho_basis->npw);

        { // it = 0
            const double dzv = static_cast<double>(GlobalC::ucell.atoms[0].ncpp.zv);
            for (int ig = igb; ig < ig_end; ++ig)
            { // initialize aux
                aux[ig] = dzv * conj(GlobalC::sf.strucFac(0, ig));
            }
        }
        for (int it = 1; it < GlobalC::ucell.ntype; it++)
        {
            const double dzv = static_cast<double>(GlobalC::ucell.atoms[it].ncpp.zv);
            for (int ig = igb; ig < ig_end; ++ig)
            { // accumulate aux
                aux[ig] += dzv * conj(GlobalC::sf.strucFac(it, ig));
            }
        }
    }

    // calculate total ionic charge
    double charge = 0.0;
    for (int it = 0; it < GlobalC::ucell.ntype; it++)
    {
        charge += GlobalC::ucell.atoms[it].na * GlobalC::ucell.atoms[it].ncpp.zv; // mohan modify 2007-11-7
    }

    double alpha = 1.1;
    double upperbound;
    do
    {
        alpha -= 0.10;
        // choose alpha in order to have convergence in the sum over G
        // upperbound is a safe upper bound for the error in the sum over G

        if (alpha <= 0.0)
        {
            ModuleBase::WARNING_QUIT("ewald", "Can't find optimal alpha.");
        }
        upperbound = 2.0 * charge * charge * sqrt(2.0 * alpha / ModuleBase::TWO_PI)
                     * erfc(sqrt(GlobalC::ucell.tpiba2 * rho_basis->ggecut / 4.0 / alpha));
    } while (upperbound > 1.0e-6);
    //	std::cout << " GlobalC::en.alpha = " << alpha << std::endl;
    //	std::cout << " upperbound = " << upperbound << std::endl;

#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for (int ig = 0; ig < rho_basis->npw; ig++)
    {
        aux[ig] *= exp(-1.0 * rho_basis->gg[ig] * GlobalC::ucell.tpiba2 / alpha / 4.0)
                   / (rho_basis->gg[ig] * GlobalC::ucell.tpiba2);
    }

    // set pos rho_basis->ig_gge0 to zero
    if (rho_basis->ig_gge0 >= 0 && rho_basis->ig_gge0 < rho_basis->npw) {
        aux[rho_basis->ig_gge0] = std::complex<double>(0.0, 0.0);
    }

#ifdef _OPENMP
#pragma omp parallel
{
    int num_threads = omp_get_num_threads();
    int thread_id = omp_get_thread_num();
#else
    int num_threads = 1;
    int thread_id = 0;
#endif

    /* Here is task distribution for multi-thread,
        0. atom will be iterated both in main nat loop and the loop in `if (rho_basis->ig_gge0 >= 0)`.
            To avoid syncing, we must calculate work range of each thread by our self
        1. Calculate the iat range [iat_beg, iat_end) by each thread
            a. when it is single thread stage, [iat_beg, iat_end) will be [0, nat)
        2. each thread iterate atoms form `iat_beg` to `iat_end-1`
    */
    int iat_beg, iat_end;
    ModuleBase::TASK_DIST_1D(num_threads, thread_id, GlobalC::ucell.nat, iat_beg, iat_end);
    iat_end = iat_beg + iat_end;

    int it_beg = (iat_beg < iat_end) ? GlobalC::ucell.iat2it[iat_beg] : GlobalC::ucell.ntype;
    int ia_beg = (iat_beg < iat_end) ? GlobalC::ucell.iat2ia[iat_beg] : 0;

    int iat = iat_beg;
    int it = it_beg;
    int ia = ia_beg;

    // preprocess ig_gap for skipping the ig point
    int ig_gap = (rho_basis->ig_gge0 >= 0 && rho_basis->ig_gge0 < rho_basis->npw) ? rho_basis->ig_gge0 : -1;

    double it_fact = 0.;
    int last_it = -1;

    // iterating atoms
    while (iat < iat_end)
    {
        if (it != last_it)
        { // calculate it_tact when it is changed
            it_fact = GlobalC::ucell.atoms[it].ncpp.zv * ModuleBase::e2 * GlobalC::ucell.tpiba
                                       * ModuleBase::TWO_PI / GlobalC::ucell.omega * fact;
            last_it = it;
        }

        const auto ig_loop = [&](int ig_beg, int ig_end)
        {
            for (int ig = ig_beg; ig < ig_end; ig++)
            {
                const ModuleBase::Vector3<double> gcar = rho_basis->gcar[ig];
                const double arg = ModuleBase::TWO_PI * (gcar * GlobalC::ucell.atoms[it].tau[ia]);
                double sumnb = -cos(arg) * aux[ig].imag() + sin(arg) * aux[ig].real();
                forceion(iat, 0) += gcar[0] * sumnb;
                forceion(iat, 1) += gcar[1] * sumnb;
                forceion(iat, 2) += gcar[2] * sumnb;
            }
        };

        // skip ig_gge0 point by separating ig loop into two part
        ig_loop(0, ig_gap);
        ig_loop(ig_gap + 1, rho_basis->npw);

        forceion(iat, 0) *= it_fact;
        forceion(iat, 1) *= it_fact;
        forceion(iat, 2) *= it_fact;

        ++iat;
        if (++ia == GlobalC::ucell.atoms[it].na)
        {
            ia = 0;
            ++it;
        }
    }

    // means that the processor contains G=0 term.
    if (rho_basis->ig_gge0 >= 0)
    {
        double rmax = 5.0 / (sqrt(alpha) * GlobalC::ucell.lat0);
        int nrm = 0;

        // output of rgen: the number of vectors in the sphere
        const int mxr = 50;
        // the maximum number of R vectors included in r
        ModuleBase::Vector3<double>* r = new ModuleBase::Vector3<double>[mxr];
        double* r2 = new double[mxr];
        ModuleBase::GlobalFunc::ZEROS(r2, mxr);
        int* irr = new int[mxr];
        ModuleBase::GlobalFunc::ZEROS(irr, mxr);
        // the square modulus of R_j-tau_s-tau_s'

        int iat1 = iat_beg;
        int T1 = it_beg;
        int I1 = ia_beg;
        const double sqa = sqrt(alpha);
        const double sq8a_2pi = sqrt(8.0 * alpha / ModuleBase::TWO_PI);

        // iterating atoms.
        // do not need to sync threads because task range of each thread is isolated
        while (iat1 < iat_end)
        {
            int iat2 = 0; // mohan fix bug 2011-06-07
            for (int T2 = 0; T2 < GlobalC::ucell.ntype; T2++)
            {
                for (int I2 = 0; I2 < GlobalC::ucell.atoms[T2].na; I2++)
                {
                    if (iat1 != iat2)
                    {
                        ModuleBase::Vector3<double> d_tau
                            = GlobalC::ucell.atoms[T1].tau[I1] - GlobalC::ucell.atoms[T2].tau[I2];
                        H_Ewald_pw::rgen(d_tau, rmax, irr, GlobalC::ucell.latvec, GlobalC::ucell.G, r, r2, nrm);

                        for (int n = 0; n < nrm; n++)
                        {
                            const double rr = sqrt(r2[n]) * GlobalC::ucell.lat0;

                            double factor
                                = GlobalC::ucell.atoms[T1].ncpp.zv * GlobalC::ucell.atoms[T2].ncpp.zv * ModuleBase::e2
                                  / (rr * rr)
                                  * (erfc(sqa * rr) / rr + sq8a_2pi * exp(-alpha * rr * rr))
                                  * GlobalC::ucell.lat0;

                            forceion(iat1, 0) -= factor * r[n].x;
                            forceion(iat1, 1) -= factor * r[n].y;
                            forceion(iat1, 2) -= factor * r[n].z;
                        }
                    }

                    ++iat2;
                }
            } // atom b

            ++iat1;
            if (++I1 == GlobalC::ucell.atoms[T1].na) {
                I1 = 0;
                ++T1;
            }
        }

        delete[] r;
        delete[] r2;
        delete[] irr;
    }
#ifdef _OPENMP
}
#endif

    Parallel_Reduce::reduce_double_pool(forceion.c, forceion.nr * forceion.nc);

    // this->print(GlobalV::ofs_running, "ewald forces", forceion);

    ModuleBase::timer::tick("Forces", "cal_force_ew");

    delete[] aux;

    return;
}

void Forces::cal_force_cc(ModuleBase::matrix& forcecc, ModulePW::PW_Basis* rho_basis, const Charge* const chr)
{
    ModuleBase::TITLE("Forces", "cal_force_cc");
    // recalculate the exchange-correlation potential.
    ModuleBase::TITLE("Forces", "cal_force_cc");
    ModuleBase::timer::tick("Forces", "cal_force_cc");

    int total_works = 0;
    // cal total works for skipping preprocess
    for (int it = 0; it < GlobalC::ucell.ntype; ++it)
    {
        if (GlobalC::ucell.atoms[it].ncpp.nlcc)
        {
            total_works += GlobalC::ucell.atoms[it].na;
        }
    }
    if (total_works == 0)
    {
        ModuleBase::timer::tick("Forces", "cal_force_cc");
        return;
    }

    ModuleBase::matrix v(GlobalV::NSPIN, rho_basis->nrxx);

    if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
    {
#ifdef USE_LIBXC
        const auto etxc_vtxc_v = XC_Functional::v_xc_meta(rho_basis->nrxx,
                                                          rho_basis->nxyz,
                                                          GlobalC::ucell.omega,
                                                          chr);

        GlobalC::en.etxc = std::get<0>(etxc_vtxc_v);
        GlobalC::en.vtxc = std::get<1>(etxc_vtxc_v);
        v = std::get<2>(etxc_vtxc_v);
#else
        ModuleBase::WARNING_QUIT("cal_force_cc", "to use mGGA, compile with LIBXC");
#endif
    }
    else
    {
        const auto etxc_vtxc_v = XC_Functional::v_xc(rho_basis->nrxx,
                                                     rho_basis->nxyz,
                                                     GlobalC::ucell.omega,
                                                     chr);

        GlobalC::en.etxc = std::get<0>(etxc_vtxc_v);
        GlobalC::en.vtxc = std::get<1>(etxc_vtxc_v);
        v = std::get<2>(etxc_vtxc_v);
    }

    const ModuleBase::matrix vxc = v;
    std::complex<double>* psiv = new std::complex<double>[rho_basis->nmaxgr];
    if (GlobalV::NSPIN == 1 || GlobalV::NSPIN == 4)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for (int ir = 0; ir < rho_basis->nrxx; ir++)
        {
            psiv[ir] = std::complex<double>(vxc(0, ir), 0.0);
        }
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for (int ir = 0; ir < rho_basis->nrxx; ir++)
        {
            psiv[ir] = 0.5 * (vxc(0, ir) + vxc(1, ir));
        }
    }

    // to G space
    rho_basis->real2recip(psiv, psiv);

#ifdef _OPENMP
#pragma omp parallel
{
    int num_threads = omp_get_num_threads();
    int thread_id = omp_get_thread_num();
#else
    int num_threads = 1;
    int thread_id = 0;
#endif

    // psiv contains now Vxc(G)
    double* rhocg = new double[rho_basis->ngg];
    ModuleBase::GlobalFunc::ZEROS(rhocg, rho_basis->ngg);

    /* Here is task distribution for multi-thread,
        0. Consider for load balancing, we distribute `total_works` instead of `nat`.
            So simply use `#pragma omp parallel for` is not suitable
        1. Calculate the work range [work_beg, work_end) by each thread
            a. when it is single thread stage, [work_beg, work_end) will be [0, total_works)
        2. each thread iterate atoms form `work_beg` to `work_end-1`
    */
    int work, work_end;
    ModuleBase::TASK_DIST_1D(num_threads, thread_id, total_works, work, work_end);
    work_end = work + work_end;

    int it = 0;
    int ia = 0;
    // We have to map the work beginning position to `it` and `ia` beginning position of this thread
    for (int work_off = 0; it < GlobalC::ucell.ntype; it++)
    {
        if (GlobalC::ucell.atoms[it].ncpp.nlcc)
        {
            if (work_off + GlobalC::ucell.atoms[it].na > work)
            {
                ia = work - work_off;
                break;
            }
            work_off += GlobalC::ucell.atoms[it].na;
        }
    }

    int last_it = -1;
    while (work < work_end)
    {
        if (it != last_it)
        {
            // call drhoc when `it` is changed
            chr->non_linear_core_correction(GlobalC::ppcell.numeric,
                                                    GlobalC::ucell.atoms[it].ncpp.msh,
                                                    GlobalC::ucell.atoms[it].ncpp.r,
                                                    GlobalC::ucell.atoms[it].ncpp.rab,
                                                    GlobalC::ucell.atoms[it].ncpp.rho_atc,
                                                    rhocg,
                                                    rho_basis);
            last_it = it;
        }

        // get iat form table
        int iat = GlobalC::ucell.itia2iat(it, ia);
        for (int ig = 0; ig < rho_basis->npw; ig++)
        {
            const ModuleBase::Vector3<double> gv = rho_basis->gcar[ig];
            const ModuleBase::Vector3<double> pos = GlobalC::ucell.atoms[it].tau[ia];
            const double rhocgigg = rhocg[rho_basis->ig2igg[ig]];
            const std::complex<double> psiv_conj = conj(psiv[ig]);

            const double arg = ModuleBase::TWO_PI * (gv.x * pos.x + gv.y * pos.y + gv.z * pos.z);
            const std::complex<double> expiarg = std::complex<double>(sin(arg), cos(arg));

            auto ipol0 = GlobalC::ucell.tpiba * GlobalC::ucell.omega * rhocgigg * gv.x * psiv_conj * expiarg;
            forcecc(iat, 0) += ipol0.real();

            auto ipol1 = GlobalC::ucell.tpiba * GlobalC::ucell.omega * rhocgigg * gv.y * psiv_conj * expiarg;
            forcecc(iat, 1) += ipol1.real();

            auto ipol2 = GlobalC::ucell.tpiba * GlobalC::ucell.omega * rhocgigg * gv.z * psiv_conj * expiarg;
            forcecc(iat, 2) += ipol2.real();
        }

        ++work;
        if (++ia == GlobalC::ucell.atoms[it].na)
        {
            ia = 0;
            do
            {   // search for next effective `it`
                ++it;
            }
            while (it < GlobalC::ucell.ntype && !GlobalC::ucell.atoms[it].ncpp.nlcc);
        }
    }

    delete[] rhocg;
#ifdef _OPENMP
} // omp parallel
#endif
    
    delete[] psiv; // mohan fix bug 2012-03-22
    Parallel_Reduce::reduce_double_pool(forcecc.c, forcecc.nr * forcecc.nc); // qianrui fix a bug for kpar > 1
    ModuleBase::timer::tick("Forces", "cal_force_cc");
    return;
}

#include "../module_base/complexarray.h"
#include "../module_base/complexmatrix.h"
void Forces::cal_force_nl(ModuleBase::matrix& forcenl, const ModuleBase::matrix& wg, const psi::Psi<complex<double>>* psi_in)
{
    ModuleBase::TITLE("Forces", "cal_force_nl");
    ModuleBase::timer::tick("Forces", "cal_force_nl");

    const int nkb = GlobalC::ppcell.nkb;
    if (nkb == 0)
        return; // mohan add 2010-07-25

    // dbecp: conj( -iG * <Beta(nkb,npw)|psi(nbnd,npw)> )
    ModuleBase::ComplexArray dbecp(3, GlobalV::NBANDS, nkb);
    ModuleBase::ComplexMatrix becp(GlobalV::NBANDS, nkb);

    // vkb1: |Beta(nkb,npw)><Beta(nkb,npw)|psi(nbnd,npw)>
    ModuleBase::ComplexMatrix vkb1(nkb, GlobalC::wf.npwx);

    for (int ik = 0; ik < GlobalC::kv.nks; ik++)
    {
        if (GlobalV::NSPIN == 2)
            GlobalV::CURRENT_SPIN = GlobalC::kv.isk[ik];
        const int nbasis = GlobalC::kv.ngk[ik];
        // generate vkb
        if (GlobalC::ppcell.nkb > 0)
        {
            GlobalC::ppcell.getvnl(ik, GlobalC::ppcell.vkb);
        }

        // get becp according to wave functions and vkb
        // important here ! becp must set zero!!
        // vkb: Beta(nkb,npw)
        // becp(nkb,nbnd): <Beta(nkb,npw)|psi(nbnd,npw)>
        becp.zero_out();
        psi_in[0].fix_k(ik);
        char transa = 'C';
        char transb = 'N';
        ///
        /// only occupied band should be calculated.
        ///
        int nbands_occ = GlobalV::NBANDS;
        while (wg(ik, nbands_occ - 1) < ModuleBase::threshold_wg)
        {
            nbands_occ--;
        }
        int npm = GlobalV::NPOL * nbands_occ;
        zgemm_(&transa,
               &transb,
               &nkb,
               &npm,
               &nbasis,
               &ModuleBase::ONE,
               GlobalC::ppcell.vkb.c,
               &GlobalC::wf.npwx,
               psi_in[0].get_pointer(),
               &GlobalC::wf.npwx,
               &ModuleBase::ZERO,
               becp.c,
               &nkb);
        Parallel_Reduce::reduce_complex_double_pool(becp.c, becp.size);

        // out.printcm_real("becp",becp,1.0e-4);
        //  Calculate the derivative of beta,
        //  |dbeta> =  -ig * |beta>
        dbecp.zero_out();
        for (int ipol = 0; ipol < 3; ipol++)
        {
            for (int i = 0; i < nkb; i++)
            {
                std::complex<double>* pvkb1 = &vkb1(i, 0);
                std::complex<double>* pvkb = &GlobalC::ppcell.vkb(i, 0);
                if (ipol == 0)
                {
                    for (int ig = 0; ig < nbasis; ig++)
                        pvkb1[ig] = pvkb[ig] * ModuleBase::NEG_IMAG_UNIT * GlobalC::wfcpw->getgcar(ik, ig)[0];
                }
                if (ipol == 1)
                {
                    for (int ig = 0; ig < nbasis; ig++)
                        pvkb1[ig] = pvkb[ig] * ModuleBase::NEG_IMAG_UNIT * GlobalC::wfcpw->getgcar(ik, ig)[1];
                }
                if (ipol == 2)
                {
                    for (int ig = 0; ig < nbasis; ig++)
                        pvkb1[ig] = pvkb[ig] * ModuleBase::NEG_IMAG_UNIT * GlobalC::wfcpw->getgcar(ik, ig)[2];
                }
            }
            std::complex<double>* pdbecp = &dbecp(ipol, 0, 0);
            zgemm_(&transa,
                &transb,
                &nkb,
                &npm,
                &nbasis,
                &ModuleBase::ONE,
                vkb1.c,
                &GlobalC::wf.npwx,
                psi_in[0].get_pointer(),
                &GlobalC::wf.npwx,
                &ModuleBase::ZERO,
                pdbecp,
                &nkb);
        }// end ipol

//		don't need to reduce here, keep dbecp different in each processor,
//		and at last sum up all the forces.
//		Parallel_Reduce::reduce_complex_double_pool( dbecp.ptr, dbecp.ndata);

//		double *cf = new double[GlobalC::ucell.nat*3];
//		ModuleBase::GlobalFunc::ZEROS(cf, GlobalC::ucell.nat);
		for (int ib=0; ib<nbands_occ; ib++)
		{
			double fac = wg(ik, ib) * 2.0 * GlobalC::ucell.tpiba;
        	int iat = 0;
        	int sum = 0;
			for (int it=0; it<GlobalC::ucell.ntype; it++)
			{
				const int Nprojs = GlobalC::ucell.atoms[it].ncpp.nh;
				for (int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
				{
					for (int ip=0; ip<Nprojs; ip++)
					{
						double ps = GlobalC::ppcell.deeq(GlobalV::CURRENT_SPIN, iat, ip, ip) ;
						const int inkb = sum + ip; 
						//out<<"\n ps = "<<ps;

						for (int ipol=0; ipol<3; ipol++)
						{
							const double dbb = ( conj( dbecp( ipol, ib, inkb) ) * becp( ib, inkb) ).real();
							forcenl(iat, ipol) = forcenl(iat, ipol) - ps * fac * dbb;
							//cf[iat*3+ipol] += ps * fac * dbb;
						}
					}

                    if(GlobalC::ppcell.multi_proj)
                    {
                        for (int ip=0; ip<Nprojs; ip++)
                        {
                            const int inkb = sum + ip;
                            //for (int ip2=0; ip2<Nprojs; ip2++)
                            for (int ip2=0; ip2<Nprojs; ip2++)
                            {
                                if ( ip != ip2 )
                                {
                                    const int jnkb = sum + ip2;
                                    double ps = GlobalC::ppcell.deeq(GlobalV::CURRENT_SPIN, iat, ip, ip2) ;

                                    for (int ipol=0; ipol<3; ipol++)
                                    {
                                        const double dbb = ( conj( dbecp( ipol, ib, inkb) ) * becp( ib, jnkb) ).real();
                                        forcenl(iat, ipol) = forcenl(iat, ipol) - ps * fac * dbb;
                                    }
                                }
                            }
                        }
                    }

					++iat;
					sum+=Nprojs;
				}
			} //end it
		} //end band
    }// end ik

    // sum up forcenl from all processors
    Parallel_Reduce::reduce_double_all(forcenl.c, forcenl.nr * forcenl.nc);
    //  this->print(GlobalV::ofs_running, "nonlocal forces", forcenl);
    ModuleBase::timer::tick("Forces", "cal_force_nl");
    return;
}

void Forces::cal_force_scc(ModuleBase::matrix& forcescc, ModulePW::PW_Basis* rho_basis)
{
    ModuleBase::TITLE("Forces", "cal_force_scc");
    ModuleBase::timer::tick("Forces", "cal_force_scc");

    //for orbital free case
    if(!GlobalC::en.vnew_exist)
    {
        return;
    }

    std::complex<double>* psic = new std::complex<double>[rho_basis->nmaxgr];

    ModuleBase::matrix& v_current = GlobalC::en.vnew;
    const int nrxx = v_current.nc;
    const int nspin = v_current.nr;

    if (nspin == 1 || nspin == 4)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for (int ir = 0; ir < nrxx; ir++)
        {
            psic[ir] = v_current(0, ir);
        }
    }
    else
    {
        int isup = 0;
        int isdw = 1;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for (int ir = 0; ir < nrxx; ir++)
        {
            psic[ir] = (v_current(isup, ir) + v_current(isdw, ir)) * 0.5;
        }
    }
    //delete vnew memory
    v_current.create(0,0);
    GlobalC::en.vnew_exist = false;

    int ndm = 0;

    for (int it = 0; it < GlobalC::ucell.ntype; it++)
    {
        if (ndm < GlobalC::ucell.atoms[it].ncpp.msh)
        {
            ndm = GlobalC::ucell.atoms[it].ncpp.msh;
        }
    }

    rho_basis->real2recip(psic, psic);

    int igg0 = 0;
    const int ig0 = rho_basis->ig_gge0;
    const int ig_gap = (ig0 >= 0 && ig0 < rho_basis->npw) ? ig0 : -1;
    if (rho_basis->gg_uniq[0] < 1.0e-8)
        igg0 = 1;

#ifdef _OPENMP
#pragma omp parallel
{
#endif

    // thread local work space
    double *aux = new double[ndm];
    ModuleBase::GlobalFunc::ZEROS(aux, ndm);
    double* rhocgnt = new double[rho_basis->ngg];
    ModuleBase::GlobalFunc::ZEROS(rhocgnt, rho_basis->ngg);

    double fact = 2.0;
    int last_it = -1;

#ifdef _OPENMP
// use no wait to avoid syncing
#pragma omp for nowait
#endif
    for (int iat = 0; iat < GlobalC::ucell.nat; ++iat)
    {
        int it = GlobalC::ucell.iat2it[iat];
        int ia = GlobalC::ucell.iat2ia[iat];

         // initialize rhocgnt when `it` is changed
        if (it != last_it)
        {
            //		Here we compute the G.ne.0 term
            const int mesh = GlobalC::ucell.atoms[it].ncpp.msh;

            for (int ig = igg0; ig < rho_basis->ngg; ++ig)
            {
                const double gx = sqrt(rho_basis->gg_uniq[ig]) * GlobalC::ucell.tpiba;
                for (int ir = 0; ir < mesh; ir++)
                {
                    if (GlobalC::ucell.atoms[it].ncpp.r[ir] < 1.0e-8)
                    {
                        aux[ir] = GlobalC::ucell.atoms[it].ncpp.rho_at[ir];
                    }
                    else
                    {
                        const double gxx = gx * GlobalC::ucell.atoms[it].ncpp.r[ir];
                        aux[ir] = GlobalC::ucell.atoms[it].ncpp.rho_at[ir] * sin(gxx) / gxx;
                    }
                }
                ModuleBase::Integral::Simpson_Integral(mesh, aux, GlobalC::ucell.atoms[it].ncpp.rab, rhocgnt[ig]);
            }

            // record it
            last_it = it;
        }

        const ModuleBase::Vector3<double> pos = GlobalC::ucell.atoms[it].tau[ia];

        const auto ig_loop = [&](int ig_beg, int ig_end)
        {
            for (int ig = ig_beg; ig < ig_end; ++ig)
            {
                const ModuleBase::Vector3<double> gv = rho_basis->gcar[ig];
                const double rhocgntigg = rhocgnt[GlobalC::rhopw->ig2igg[ig]];
                const double arg = ModuleBase::TWO_PI * (gv * pos);
                const std::complex<double> cpm = std::complex<double>(sin(arg), cos(arg)) * conj(psic[ig]);

                forcescc(iat, 0) += fact * rhocgntigg * GlobalC::ucell.tpiba * gv.x * cpm.real();
                forcescc(iat, 1) += fact * rhocgntigg * GlobalC::ucell.tpiba * gv.y * cpm.real();
                forcescc(iat, 2) += fact * rhocgntigg * GlobalC::ucell.tpiba * gv.z * cpm.real();
            }
        };

        ig_loop(0, ig_gap);
        ig_loop(ig_gap + 1, rho_basis->npw);

        // std::cout << " forcescc = " << forcescc(iat,0) << " " << forcescc(iat,1) << " " <<
        // forcescc(iat,2) << std::endl;
    }

    delete[] aux;
    delete[] rhocgnt;

#ifdef _OPENMP
}
#endif

    Parallel_Reduce::reduce_double_pool(forcescc.c, forcescc.nr * forcescc.nc);

    delete[] psic; // mohan fix bug 2012-03-22

    ModuleBase::timer::tick("Forces", "cal_force_scc");
    return;
}
