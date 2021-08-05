#ifndef MODEC_SOLVER_H_
#define MODEC_SOLVER_H_

#include "sparse_lu.h"
#include "spdlog.h"

namespace MODEC {
    const std::map<std::string, int> SolverIndicator{
        {"PFDCRAM", 1},
        {"IPFCRAM", 2},
        {"PFDPRAM", 3},
        {"IPFPRAM", 4}
    };
};


class ModecSolver {
public:
    // singleton pattern to construct modec solver
    static ModecSolver& getModecSolver();

    void SetSolverParams(const std::string& _name, const int& _major_order, const int& _minor_order = 0);

    // perform modec solver
    void Execute(SparseMatrix<double>& _burnup_matrix, SymbolicSpMat& _sym_spmat, double& _time, std::vector<double>& _nuclide_densities);

private:
    ModecSolver() {};
    ~ModecSolver() {};
    ModecSolver(const ModecSolver&);
    ModecSolver& operator=(const ModecSolver&);

    int solver_id_{ 2 };
    int major_order_{ 48 };
    int minor_order_{ 0 };
};

// pfd-cram solver
template<int ORDER>
class PfdCRAM {
public:
    static void Solve(SparseMatrix<double>& _burnup_matrix, SymbolicSpMat& _sym_spmat, const double& _time, std::vector<double>& _nuclide_densities) {
        std::vector<std::complex<double> > residues(ORDER / 2); // namely "alpha" in previous version
        std::vector<std::complex<double> > poles(ORDER / 2);    // namely "theta" in previous version

        double residue_inf = 0.0;

        CalcResiduesPoles(residues, poles, residue_inf);

        int nucl_dim = _nuclide_densities.size();
        std::vector<std::complex<double> > nucl_dens_temp(nucl_dim);

        std::vector<double> nucl_dens_final;
        nucl_dens_final.resize(nucl_dim);

        SparseLU<std::complex<double> > solver;

        for (int j = 0; j < ORDER / 2; ++j) {
            for (int i = 0; i < nucl_dim; ++i) {
                nucl_dens_temp[i] = _nuclide_densities[i] * residues[j];
            }

            solver.SolverLU(_burnup_matrix * _time - poles[j], _sym_spmat, nucl_dens_temp);

            for (int i = 0; i < nucl_dim; ++i) {
                nucl_dens_final[i] += 2.0 * nucl_dens_temp[i].real();
            }
        }
        for (unsigned int i = 0; i < nucl_dim; ++i) {
            _nuclide_densities[i] = nucl_dens_final[i] + _nuclide_densities[i] * residue_inf;
        }
    }
private:
    static void CalcResiduesPoles(std::vector<std::complex<double> >& _residues, std::vector<std::complex<double> >& _poles, double& _residue_inf) {
        switch (ORDER) {
        case 14:
            _residues = {
                std::complex<double>(-7.1542880635890672853E-05, +1.4361073349541300111E-04),
                std::complex<double>(+9.4390253107361688779E-03, -1.7184791958483017511E-02),
                std::complex<double>(-3.7636003878226968717E-01, +3.3518347029450104214E-01),
                std::complex<double>(-2.3498232091082701191E+01, -5.8083591297142074004E+00),
                std::complex<double>(+4.6933274488831293047E+01, +4.5643649768827760791E+01),
                std::complex<double>(-2.7875161940145646468E+01, -1.0214733999056451434E+02),
                std::complex<double>(+4.8071120988325088907E+00, -1.3209793837428723881E+00)
            };

            _poles = {
                std::complex<double>(-8.8977731864688888199E+00, +1.6630982619902085304E+01),
                std::complex<double>(-3.7032750494234480603E+00, +1.3656371871483268171E+01),
                std::complex<double>(-0.2087586382501301251E+00, +1.0991260561901260913E+01),
                std::complex<double>(+3.9933697105785685194E+00, +6.0048316422350373178E+00),
                std::complex<double>(+5.0893450605806245066E+00, +3.5888240290270065102E+00),
                std::complex<double>(+5.6231425727459771248E+00, +1.1940690463439669766E+00),
                std::complex<double>(+2.2697838292311127097E+00, +8.4617379730402214019E+00)
            };

            _residue_inf = 1.8321743782540412751E-14;

            break;

        case 16:
            _residues = {
                std::complex<double>(-5.0901521865224915650E-07, -2.4220017652852287970E-05),
                std::complex<double>(+2.1151742182466030907E-04, +4.3892969647380673918E-03),
                std::complex<double>(+1.1339775178483930527E+02, +1.0194721704215856450E+02),
                std::complex<double>(+1.5059585270023467528E+01, -5.7514052776421819979E+00),
                std::complex<double>(-6.4500878025539646595E+01, -2.2459440762652096056E+02),
                std::complex<double>(-1.4793007113557999718E+00, +1.7686588323782937906E+00),
                std::complex<double>(-6.2518392463207918892E+01, -1.1190391094283228480E+01),
                std::complex<double>(+4.1023136835410021273E-02, -1.5743466173455468191E-01)
            };

            _poles = {
                std::complex<double>(-1.0843917078696988026E+01, +1.9277446167181652284E+01),
                std::complex<double>(-5.2649713434426468895E+00, +1.6220221473167927305E+01),
                std::complex<double>(+5.9481522689511774808E+00, +3.5874573620183222829E+00),
                std::complex<double>(+3.5091036084149180974E+00, +8.4361989858843750826E+00),
                std::complex<double>(+6.4161776990994341923E+00, +1.1941223933701386874E+00),
                std::complex<double>(+1.4193758971856659786E+00, +1.0925363484496722585E+01),
                std::complex<double>(+4.9931747377179963991E+00, +5.9968817136039422260E+00),
                std::complex<double>(-1.4139284624888862114E+00, +1.3497725698892745389E+01)
            };

            _residue_inf = 2.1248537104952237488E-16;

            break;

        default:
            std::cerr << "Error: the order of PFD-CRAM method must be either 14 or 16 !! \n";

            // RuntimeInformation::LogError("Position: void PfdCramSolver::CalcResiduesPoles; \n Error: the order of the CRAM method must be either 14 or 16 !!");

            break;
        }
    }
};

// ipf-cram solver
template<int ORDER>
class IpfCRAM {
public:
    static void Solve(SparseMatrix<double>& _burnup_matrix, SymbolicSpMat& _sym_spmat, const double& _time, std::vector<double>& _nuclide_densities) {
        std::vector<std::complex<double> > residues(ORDER / 2); // namely "alpha" in previous version
        std::vector<std::complex<double> > poles(ORDER / 2);    // namely "theta" in previous version

        double residue_inf = 0.0;

        CalcResiduesPoles(residues, poles, residue_inf);

        int nucl_dim = _nuclide_densities.size();

        std::vector<std::complex<double> > nucl_dens_temp(nucl_dim);

        std::vector<double> nucl_dens_final(_nuclide_densities);

        SparseLU<std::complex<double> > solver;

        for (int j = 0; j < ORDER / 2; ++j) {
            //_burnup_matrix->DoConcreteLuFactorization(poles[j], _time, nucl_dens_final, nucl_dens_temp);
            nucl_dens_temp = solver.SolverLU(_burnup_matrix * _time - poles[j], _sym_spmat, nucl_dens_final); // At - theta

            for (int i = 0; i < nucl_dim; ++i) {
                nucl_dens_final[i] += 2.0 * (residues[j] * nucl_dens_temp[i]).real();
            }
        }
        for (unsigned int i = 0; i < nucl_dim; ++i) {
            _nuclide_densities[i] = nucl_dens_final[i] * residue_inf;
        }
    }
private:
    static void CalcResiduesPoles(std::vector<std::complex<double> >& _residues, std::vector<std::complex<double> >& _poles, double& _residue_inf) {
        switch (ORDER) {
        case 32:
            _residues = {
                std::complex<double>(+2.531738087291248E+03, -2.175335564251554E+05),
                std::complex<double>(+1.124951925994460E+02, -2.293141996821969E+02),
                std::complex<double>(+1.928222545500035E+02, -5.268292754604315E+02),
                std::complex<double>(+4.159057536149641E+04, -3.509499779111824E+05),
                std::complex<double>(+4.981965659174993E+03, -5.519940772045004E+05),
                std::complex<double>(+1.498320382271818E+02, -6.815792364464349E+02),
                std::complex<double>(+6.462582753425457E+02, -6.286003508366936E+03),
                std::complex<double>(+3.108473530705140E+02, -2.227478875866931E+04),
                std::complex<double>(+2.885604705807475E+02, -1.308904072042900E+03),
                std::complex<double>(+7.047430374921731E+01, -3.234996756134937E+02),
                std::complex<double>(+7.344336306595115E+01, -2.581272076901578E+02),
                std::complex<double>(+4.519923142224831E+01, -4.759858266475396E+00),
                std::complex<double>(+5.271307338359870E+01, -2.583918268027449E+01),
                std::complex<double>(+4.983259420279240E+01, -1.542812979380035E+01),
                std::complex<double>(+6.057642473653785E+01, -6.941674585178426E+02),
                std::complex<double>(+4.534413157137819E+01, -3.114070953092643E+01)
            };

            _poles = {
                std::complex<double>(+9.093582328296485E+00, +1.321250259877626E+01),
                std::complex<double>(-1.727168130828052E+00, +2.572302624565582E+01),
                std::complex<double>(+1.166498099877330E+00, +2.314228200408637E+01),
                std::complex<double>(+5.780372991371440E+00, +1.811790953587012E+01),
                std::complex<double>(+1.201313618215791E+01, +5.977830703724750E+00),
                std::complex<double>(+7.585235004412852E+00, +1.565385209991484E+01),
                std::complex<double>(+3.652646888760184E+00, +2.061118862125454E+01),
                std::complex<double>(+1.272430588935007E+01, +1.194256446871873E+00),
                std::complex<double>(-5.098669890503975E+00, +2.837060507799692E+01),
                std::complex<double>(+1.129472112530614E+01, +8.378353071561875E+00),
                std::complex<double>(+1.032515517992622E+01, +1.078873149419210E+01),
                std::complex<double>(-2.734798356101249E+01, +4.069231229671756E+01),
                std::complex<double>(-1.377054550962031E+01, +3.399214933428828E+01),
                std::complex<double>(-1.958248015332671E+01, +3.710329224494451E+01),
                std::complex<double>(+1.248804035352348E+01, +3.584068696364335E+00),
                std::complex<double>(-9.054506283982660E+00, +3.111173004087709E+01)
            };

            _residue_inf = 6.932444346272945E-32;
            break;

        case 48:
            _residues = {
                std::complex<double>(+6.387380733878774E+02, -6.743912502859256E+02),
                std::complex<double>(+1.909896179065730E+02, -3.973203432721332E+02),
                std::complex<double>(+4.236195226571914E+02, -2.041233768918671E+03),
                std::complex<double>(+4.645770595258726E+02, -1.652917287299683E+03),
                std::complex<double>(+7.765163276752433E+02, -1.783617639907328E+04),
                std::complex<double>(+1.907115136768522E+03, -5.887068595142284E+04),
                std::complex<double>(+2.909892685603256E+03, -9.953255345514560E+03),
                std::complex<double>(+1.944772206620450E+02, -1.427131226068449E+03),
                std::complex<double>(+1.382799786972332E+05, -3.256885197214938E+06),
                std::complex<double>(+5.628442079602433E+03, -2.924284515884309E+04),
                std::complex<double>(+2.151681283794220E+02, -1.121774011188224E+03),
                std::complex<double>(+1.324720240514420E+03, -6.370088443140973E+04),
                std::complex<double>(+1.617548476343347E+04, -1.008798413156542E+06),
                std::complex<double>(+1.112729040439685E+02, -8.837109731680418E+01),
                std::complex<double>(+1.074624783191125E+02, -1.457246116408180E+02),
                std::complex<double>(+8.835727765158191E+01, -6.388286188419360E+01),
                std::complex<double>(+9.354078136054179E+01, -2.195424319460237E+02),
                std::complex<double>(+9.418142823531573E+01, -6.719055740098035E+02),
                std::complex<double>(+1.040012390717851E+02, -1.693747595553868E+02),
                std::complex<double>(+6.861882624343235E+01, -1.177598523430493E+01),
                std::complex<double>(+8.766654491283722E+01, -4.596464999363902E+03),
                std::complex<double>(+1.056007619389650E+02, -1.738294585524067E+03),
                std::complex<double>(+7.738987569039419E+01, -4.311715386228984E+01),
                std::complex<double>(+1.041366366475571E+02, -2.777743732451969E+02)
            };

            _poles = {
                std::complex<double>(-4.465731934165702E+01, +6.233225190695437E+01),
                std::complex<double>(-5.284616241568964E+00, +4.057499381311059E+01),
                std::complex<double>(-8.867715667624458E+00, +4.325515754166724E+01),
                std::complex<double>(+3.493013124279215E+00, +3.281615453173585E+01),
                std::complex<double>(+1.564102508858634E+01, +1.558061616372237E+01),
                std::complex<double>(+1.742097597385893E+01, +1.076629305714420E+01),
                std::complex<double>(-2.834466755180654E+01, +5.492841024648724E+01),
                std::complex<double>(+1.661569367939544E+01, +1.316994930024688E+01),
                std::complex<double>(+8.011836167974721E+00, +2.780232111309410E+01),
                std::complex<double>(-2.056267541998229E+00, +3.794824788914354E+01),
                std::complex<double>(+1.449208170441839E+01, +1.799988210051809E+01),
                std::complex<double>(+1.853807176907916E+01, +5.974332563100539E+00),
                std::complex<double>(+9.932562704505182E+00, +2.532823409972962E+01),
                std::complex<double>(-2.244223871767187E+01, +5.179633600312162E+01),
                std::complex<double>(+8.590014121680897E-01, +3.536456194294350E+01),
                std::complex<double>(-1.286192925744479E+01, +4.600304902833652E+01),
                std::complex<double>(+1.164596909542055E+01, +2.287153304140217E+01),
                std::complex<double>(+1.806076684783089E+01, +8.368200580099821E+00),
                std::complex<double>(+5.870672154659249E+00, +3.029700159040121E+01),
                std::complex<double>(-3.542938819659747E+01, +5.834381701800013E+01),
                std::complex<double>(+1.901323489060250E+01, +1.194282058271408E+00),
                std::complex<double>(+1.885508331552577E+01, +3.583428564427879E+00),
                std::complex<double>(-1.734689708174982E+01, +4.883941101108207E+01),
                std::complex<double>(+1.316284237125190E+01, +2.042951874827759E+01)
            };

            _residue_inf = 2.258038182743983E-47;
            break;

        default:
            std::cerr << "Error: the order of IPF-CRAM method must be either 32 or 48 !! \n";

            //// RuntimeInformation::LogError("Position: void IpfCramSolver::CalcResiduesPoles; \n Error: the order of the IPF-CRAM method must be either 32 or 48 !!");

            break;

        }
    }
};

// pfd-pram solver
template<int MAJOR_ORDER, int MINOR_ORDER>
class PfdPRAM {
public:
    static void Solve(SparseMatrix<double>& _burnup_matrix, SymbolicSpMat& _sym_spmat, const double& _time, std::vector<double>& _nuclide_densities) {
        std::vector<std::complex<double> > residues(MAJOR_ORDER / 2); // namely "alpha" in previous version
        std::vector<std::complex<double> > poles(MAJOR_ORDER / 2);    // namely "theta" in previous version

        double residue_inf = 0.0;

        CalcResiduesPoles(residues, poles, residue_inf);

        int nucl_dim = _nuclide_densities.size();
        std::vector<std::complex<double> > nucl_dens_temp(nucl_dim);

        std::vector<double> nucl_dens_final;
        nucl_dens_final.resize(nucl_dim);

        SparseLU<std::complex<double> > solver;

        for (int j = MAJOR_ORDER / 2 - 1; j >= 0; --j) {
            for (int i = 0; i < nucl_dim; ++i) {
                nucl_dens_temp[i] = _nuclide_densities[i];
            }
            //_burnup_matrix->DoConcreteLuFactorization(poles[j], _time, nucl_dens_temp);
            solver.SolverLU(_burnup_matrix * _time - poles[j], _sym_spmat, nucl_dens_temp);

            for (int i = 0; i < nucl_dim; ++i) {
                nucl_dens_final[i] += (residues[j] * nucl_dens_temp[i]).real();
            }
        }
        _nuclide_densities = nucl_dens_final;
    }
private:
    static void CalcResiduesPoles(std::vector<std::complex<double> >& _residues, std::vector<std::complex<double> >& _poles, double& _residue_inf) {
        switch (MAJOR_ORDER) {
        case 16: {
            switch (MINOR_ORDER) {
            case 4: {
                _residues = {
                    std::complex<double>(-4105.833110178047871872696357517, -9978.7456818627129650639266575943),
                    std::complex<double>(7074.827941091885188510940145947, 2950.9995256349586360607463967575),
                    std::complex<double>(-3503.4143385356074746948053930378, 1498.5659416629983794137043284732),
                    std::complex<double>(430.52713795461152852748871243797, -1206.7333811543962924417913952527),
                    std::complex<double>(135.6479051361795847000566708487, 238.02432302877453635803050954583),
                    std::complex<double>(-32.995666254016041225776176931376, -4.6641750946273966275776165214256),
                    std::complex<double>(1.2312557109715693380621885604132, -1.3787952522362761856349073156036),
                    std::complex<double>(0.008875074023516716730209692049139, 0.025194171556626594155667678527108)
                };

                _poles = {
                    std::complex<double>(9.64137579910935864398505377985, 1.0981010822149490368561033740966),
                    std::complex<double>(9.2763242514386015896566210825134, 3.2895827751721238126596824280871),
                    std::complex<double>(8.52965610659565667298635229736, 5.466576925530935315222818508057),
                    std::complex<double>(7.3643007632255104919339755394608, 7.6182859115883767373836600961492),
                    std::complex<double>(5.7123533143166997272494675459542, 9.732015526205308892922069911635),
                    std::complex<double>(3.4501613245318985927500480131446, 11.79205726063087887084147418132),
                    std::complex<double>(0.32890188837385537793559311481635, 13.778597759288537022581692627942),
                    std::complex<double>(-4.3030734475915810964971113730993, 15.671313533282350685478368920649)
                };
            }
                  break;

            case 8: {
                _residues = {
                    std::complex<double>(-123885.60732478637914000570539081, -788788.11018923403910678542648302),
                    std::complex<double>(251575.12786960887559715888782385, 499201.45540701365270921333153376),
                    std::complex<double>(-188986.73851511562160616277182519, -192647.57083802195758899087547096),
                    std::complex<double>(76596.475788147093454106039648335, 40935.082868012145510197352457381),
                    std::complex<double>(-17184.033656053279613316653078058, -3435.641601798097967617544337095),
                    std::complex<double>(1981.3799168100933120597802610125, -143.32568478748701800223227201584),
                    std::complex<double>(-97.91430106013304928326126468782, 28.62389985039226876243858128279),
                    std::complex<double>(1.3102224493510454436838255499764, -0.44622624894776383870012486618064)
                };

                _poles = {
                    std::complex<double>(13.742665875890958130639663887724, 1.345720189669652107588108489853),
                    std::complex<double>(13.366449189282172952275768073607, 4.0375473255429475397229728068698),
                    std::complex<double>(12.597572537563548371990484765739, 6.7308791533334698280155244480057),
                    std::complex<double>(11.39935605569665519559656186862, 9.4280597347684155361864484561328),
                    std::complex<double>(9.7049674785523400464775935560337, 12.134315480013574658887761298754),
                    std::complex<double>(7.3935898086823968129855501883908, 14.861728331100349273241628075976),
                    std::complex<double>(4.2242166799829803390633153296833, 17.64129629009021381437844350719),
                    std::complex<double>(-0.42881762565105184902893766979765, 20.576995043965643935544486997952)
                };
            }
                  break;

            default: {
                std::cerr << "Error: the minor order of 16th-order PFD-PRAM method must be either 4 or 8 !! \n";
                // RuntimeInformation::LogError("Position: void PfdCramSolver::CalcResiduesPoles; \n Error: the minor order of 16th-order PFD-PRAM method must be either 4 or 8 !!");
            }
                   break;
            }
        }
               break;

        case 32: {
            switch (MINOR_ORDER) {
            case 4: {
                _residues = {
                    std::complex<double>(-26207.786112846926327164233644596, 2013.9574304715661340691232390949),
                    std::complex<double>(915096.0667448110096490157543599, -196287.55745638523689858012537801),
                    std::complex<double>(-21.285490264930621207289319377655, -16.128180421061280628326812811588),
                    std::complex<double>(149056.96205567608422522877983659, -124724.43645329372340730026159601),
                    std::complex<double>(-621029.76325101138572534436023166, -920984.0653441356105954464356361),
                    std::complex<double>(-122.31075798794265430022255005811, 196.13058092869503485984007309494),
                    std::complex<double>(0.089333451435590367291273201140177, 0.054291831296643836055117131395126),
                    std::complex<double>(-0.000000046206222308423851829062801135933, -0.000000032915332440826664936495918937387),
                    std::complex<double>(28809.503252551828407356632263684, 73755.73893332862606667624612787),
                    std::complex<double>(-340051.88914745686368198400859574, -198014.10634347516457484817409604),
                    std::complex<double>(1272.4737367909591418777172550311, -6875.9629826443059481700226196401),
                    std::complex<double>(-0.000028470676852187418514815250944518, 0.0000047922906188279753373141512382154),
                    std::complex<double>(-0.0007958277553947126811631439674361, 0.002641841278463188306731320653845),
                    std::complex<double>(-108147.41979115556466296400805495, 654228.38175172260314397969481991),
                    std::complex<double>(1.284518043902611191813520883228, -1.6802261499502779830187278362789),
                    std::complex<double>(1344.0757337430325171346574164272, 565.06006111596823848009542168512)
                };

                _poles = {
                    std::complex<double>(10.552083182602945373543744559998, 11.866991747737281587031497654451),
                    std::complex<double>(14.266047827278894324113091545026, 2.7870345574994553837123977389392),
                    std::complex<double>(3.3745062813808144164020560875694, 18.486240713334667123478513543419),
                    std::complex<double>(12.63245941275774879766481226944, 8.2976664026921784824909969156421),
                    std::complex<double>(14.444057281439322039924063677589, 0.92977355413925930729762487969671),
                    std::complex<double>(5.624250033509810558180771581522, 16.922639694208363086408355698398),
                    std::complex<double>(-2.4154155668738209572142789352004, 21.32730931662689444776021044374),
                    std::complex<double>(-17.648705066703518986028555196597, 24.233470743717409996589999806548),
                    std::complex<double>(11.698918043504880882121371733434, 10.096798762759585759134538927882),
                    std::complex<double>(13.365662333613262866767456461496, 6.4759718079259341780105301114975),
                    std::complex<double>(9.1744903734590893786139991115289, 13.6008687027368391346929042319),
                    std::complex<double>(-11.022924484632429921811464008774, 23.553226330713849754186821285781),
                    std::complex<double>(-6.2240763146972220688029619283376, 22.543994866250547125665688402376),
                    std::complex<double>(13.907984902718445211471496290129, 4.6373753934816370183417480442461),
                    std::complex<double>(0.72806817115414645484293108265729, 19.962650509055565451223081169109),
                    std::complex<double>(7.5425935894876316302114656685189, 15.28966130919649377951409825212)
                };
            }
                  break;

            case 8: {
                _residues = {
                    std::complex<double>(11307275.165746824292014188542422, 9237594.0294776047546183909315407),
                    std::complex<double>(2.5106103894463108911836597189629, -4.2978134131544908988284220422573),
                    std::complex<double>(647182.81558960788928317473228702, -30341818.391197467408941900977944),
                    std::complex<double>(607.86297944800133608312124900861, 1348.5752279049379850362734368771),
                    std::complex<double>(0.000000040760380603243897574833420632252, -0.0000023643121495258992063147768284058),
                    std::complex<double>(0.072676736986825114295411161913688, 0.097858980534346433409471444724523),
                    std::complex<double>(718108.74923731083334199622943745, -1688950.8395982116300002187403931),
                    std::complex<double>(72433545.472834132456148547684777, 17721192.247926872685099288273237),
                    std::complex<double>(-5597703.9910874664338125765550998, 1335766.7901388949689550922884038),
                    std::complex<double>(-91874.848854385860025992990707086, -4009.6385960066022781736718811269),
                    std::complex<double>(243928.95084602926116798104927943, 396734.64964469047364395357997597),
                    std::complex<double>(-0.0011526395699856839352830254615718, 0.00034444655915066609920813035303912),
                    std::complex<double>(7869.6688985175575670897399749224, -11224.501720052710962334579058846),
                    std::complex<double>(-38199806.89111815845944761290347, -80472994.382665905267436722751196),
                    std::complex<double>(-41469027.048883451765650024757009, 31625818.858549457843279942954721),
                    std::complex<double>(-108.48832293539545377867875416915, -8.1152456803955282280178573946514)
                };

                _poles = {
                    std::complex<double>(16.804372519568827486391980973765, 9.7292083066289292558950568191997),
                    std::complex<double>(1.3661602413543764359708144178536, 26.161074809918883557542844100248),
                    std::complex<double>(17.560852826236411757764433931459, 7.5794766655203362015666896734585),
                    std::complex<double>(7.2811307966066359217653249014532, 22.243674462428705013712928467104),
                    std::complex<double>(-13.927676223891141936571219058009, 31.648395252811623813856865562654),
                    std::complex<double>(-2.4995310888206177836637506887431, 28.044712824816047873410456251928),
                    std::complex<double>(14.65938628923407707267464572919, 13.988593588933133231252958058086),
                    std::complex<double>(18.490141199340384891975196574175, 3.2548348766449188728805588482706),
                    std::complex<double>(15.84156092795897839863669011429, 11.86665176396150252695880368106),
                    std::complex<double>(11.560952411847126281192633355946, 18.171405962478540462375165666209),
                    std::complex<double>(13.240346369411666492588074610952, 16.091496275461808637008819303807),
                    std::complex<double>(-7.3339837729370368495722903332941, 29.870416314615500583743774181372),
                    std::complex<double>(9.5893163866070215505355980671503, 20.223842625789631404428838786093),
                    std::complex<double>(18.673903714987714948137831741339, 1.0853727152211989477412078262598),
                    std::complex<double>(18.120544172173429490704087812923, 5.4204235869878144503289374616496),
                    std::complex<double>(4.5725232303221458414699478495491, 24.224994315300233600329341409826)
                };
            }
                  break;

            default: {
                std::cerr << "Error: the minor order of 32th-order PFD-PRAM method must be either 4 or 8 !! \n";
                // RuntimeInformation::LogError("Position: void PfdCramSolver::CalcResiduesPoles; \n Error: the minor order of 16th-order PFD-PRAM method must be either 4 or 8 !!");
            }
                   break;
            }
        }
               break;

        default: {
            std::cerr << "Error: the order of PFD-PRAM method must be either 16 or 32 !! \n";
            // RuntimeInformation::LogError("Position: void PfdCramSolver::CalcResiduesPoles; \n Error: the order of the CRAM method must be either 14 or 16 !!");
        }
               break;

        }
    }
};

// ipf-pram solver
template<int MAJOR_ORDER, int MINOR_ORDER>
class IpfPRAM {
public:
    static void Solve(SparseMatrix<double>& _burnup_matrix, SymbolicSpMat& _sym_spmat, const double& _time, std::vector<double>& _nuclide_densities) {
        std::vector<std::complex<double> > residues(MAJOR_ORDER / 2); // namely "alpha" in previous version
        std::vector<std::complex<double> > poles(MAJOR_ORDER / 2);    // namely "theta" in previous version

        double residue_inf = 0.0;

        CalcResiduesPoles(residues, poles, residue_inf);

        int nucl_dim = _nuclide_densities.size();

        std::vector<std::complex<double> > nucl_dens_temp(nucl_dim);

        std::vector<double> nucl_dens_final(_nuclide_densities);

        SparseLU<std::complex<double> > solver;

        for (int j = MAJOR_ORDER / 2 - 1; j > MINOR_ORDER / 2 - 1; --j) {
            nucl_dens_temp = solver.SolverLU(_burnup_matrix * _time - poles[j], _sym_spmat, nucl_dens_final);
            for (unsigned int i = 0; i < nucl_dim; ++i) {
                nucl_dens_final[i] = (residues[j] * nucl_dens_temp[i]).real();
            }
        }
        for (int j = MINOR_ORDER / 2 - 1; j >= 0; --j) {
            nucl_dens_temp = solver.SolverLU(_burnup_matrix * _time - poles[j], _sym_spmat, nucl_dens_final);
            for (unsigned int i = 0; i < nucl_dim; ++i) {
                nucl_dens_final[i] += (residues[j] * nucl_dens_temp[i]).real();
            }
        }

        for (unsigned int i = 0; i < nucl_dim; ++i) {
            _nuclide_densities[i] = nucl_dens_final[i] * residue_inf;
        }
    }
private:
    static void CalcResiduesPoles(std::vector<std::complex<double> >& _residues, std::vector<std::complex<double> >& _poles, double& _residue_inf) {
        switch (MAJOR_ORDER) {
        case 16: {
            switch (MINOR_ORDER) {
            case 4: {
                _residues = {
                    std::complex<double>(27.011452627867510067844280741766, 3.4182754602096200698125532350975),
                    std::complex<double>(33.040204253697038495032682741668, -12.828284645081451352373606867956),
                    std::complex<double>(0, -8.4802844651934778437621192435461),
                    std::complex<double>(0, -10.275363795991787902401812115365),
                    std::complex<double>(0, -13.126312291310483372085626900787),
                    std::complex<double>(0, -18.292983225565349413750128454599),
                    std::complex<double>(0, -30.39899185840295799223992406689),
                    std::complex<double>(0, -91.066297647474121176236310398002)
                };

                _poles = {
                    std::complex<double>(-4.3030734475915810964971113730993, 15.671313533282350685478368920649),
                    std::complex<double>(0.32890188837385537793559311481635, 13.778597759288537022581692627942),
                    std::complex<double>(3.4501613245318985927500480131446, 11.79205726063087887084147418132),
                    std::complex<double>(5.7123533143166997272494675459542, 9.732015526205308892922069911635),
                    std::complex<double>(7.3643007632255104919339755394608, 7.6182859115883767373836600961492),
                    std::complex<double>(8.52965610659565667298635229736, 5.466576925530935315222818508057),
                    std::complex<double>(9.2763242514386015896566210825134, 3.2895827751721238126596824280871),
                    std::complex<double>(9.64137579910935864398505377985, 1.0981010822149490368561033740966)
                };

                _residue_inf = 0.871782912;
            }
                  break;

            case 8: {
                _residues = {
                    std::complex<double>(37.546000850596167417255130438190536836184225628934, 3.1904280941673704827161545935871537346811762864613),
                    std::complex<double>(45.294669151967774654993954077183228048131567039427, -14.200350156361281413143672175108026380085532677745),
                    std::complex<double>(48.242508193494493928797351424511793308523228600983, -33.767273861793252666919838512538677629982168017435),
                    std::complex<double>(46.704734487074894697948606868734489940544672538609, -57.460082844049698979106766526643652930715352533834),
                    std::complex<double>(0, -10.606636234094281771660636579425089208926706581898),
                    std::complex<double>(0, -14.856900223869711190089221431009158567805003242773),
                    std::complex<double>(0, -24.767511545280163506549068330123241590001839682805),
                    std::complex<double>(0, -74.309652755189797051002449927977595280242503386302)
                };

                _poles = {
                    std::complex<double>(-0.42881762565105184902893766979765367646300577372985, 20.576995043965643935544486997951685964138066863983),
                    std::complex<double>(4.2242166799829803390633153296832494123035621678448, 17.641296290090213814378443507189552808355058489111),
                    std::complex<double>(7.3935898086823968129855501883907575181425225452071, 14.86172833110034927324162807597599812933844665669),
                    std::complex<double>(9.7049674785523400464775935560336708127087679646545, 12.134315480013574658887761298754394600469220805759),
                    std::complex<double>(11.399356055696655195596561868620054682926415178797, 9.4280597347684155361864484561328085989596761865756),
                    std::complex<double>(12.597572537563548371990484765738961847346583590839, 6.7308791533334698280155244480056704411567978833091),
                    std::complex<double>(13.366449189282172952275768073607416171172923200192, 4.037547325542947539722972806869754963453831265708),
                    std::complex<double>(13.742665875890958130639663887723543231862231126195, 1.3457201896696521075881084898529634607795410729212)
                };

                _residue_inf = 5.189184;
            }
                  break;

            default: {
                std::cerr << "Error: the minor order of 16th-order IPF-PRAM method must be either 4 or 8 !! \n";
                //// RuntimeInformation::LogError("Position: void IpfPramSolver::CalcResiduesPoles; \n Error: the minor order of 16th-order IPF-PRAM method must be either 4 or 8 !!");
            }
                   break;
            }
        }
               break;

        case 32: {
            switch (MINOR_ORDER) {
            case 8: {
                _residues = {
                    std::complex<double>(42.580230704504434701484939832704, 17.010770869305818263787549484256),
                    std::complex<double>(54.194481088010281057034065783068, 2.2193147090650595972640619333658),
                    std::complex<double>(60.445023068598656852348348170542, -14.030816581246850170827847675406),
                    std::complex<double>(61.990203450297787121459754889301, -32.507421151675002363575587111037),
                    std::complex<double>(0, -4.1279679449435877121214289986333),
                    std::complex<double>(0, -4.4956601108736675760885912406396),
                    std::complex<double>(0, -4.9446587303087037438002719555782),
                    std::complex<double>(0, -5.5031515011268955401113253929096),
                    std::complex<double>(0, -6.2144624892646975147803795937424),
                    std::complex<double>(0, -7.1486814856865600537785173788012),
                    std::complex<double>(0, -8.4269768751195332805983147776214),
                    std::complex<double>(0, -10.278328600680251324997099814068),
                    std::complex<double>(0, -13.193523037666734783370149482051),
                    std::complex<double>(0, -18.448742684992084931379482374596),
                    std::complex<double>(0, -30.723524783868580811369418135851),
                    std::complex<double>(0, -92.134248998160969129098074174584)
                };

                _poles = {
                    std::complex<double>(-13.927676223891141936571219058009, 31.648395252811623813856865562654),
                    std::complex<double>(-7.3339837729370368495722903332941, 29.870416314615500583743774181372),
                    std::complex<double>(-2.4995310888206177836637506887431, 28.044712824816047873410456251928),
                    std::complex<double>(1.3661602413543764359708144178536, 26.161074809918883557542844100248),
                    std::complex<double>(4.5725232303221458414699478495491, 24.224994315300233600329341409826),
                    std::complex<double>(7.2811307966066359217653249014532, 22.243674462428705013712928467104),
                    std::complex<double>(9.5893163866070215505355980671503, 20.223842625789631404428838786093),
                    std::complex<double>(11.560952411847126281192633355946, 18.171405962478540462375165666209),
                    std::complex<double>(13.240346369411666492588074610952, 16.091496275461808637008819303807),
                    std::complex<double>(14.65938628923407707267464572919, 13.988593588933133231252958058086),
                    std::complex<double>(15.84156092795897839863669011429, 11.86665176396150252695880368106),
                    std::complex<double>(16.804372519568827486391980973765, 9.7292083066289292558950568191997),
                    std::complex<double>(17.560852826236411757764433931459, 7.5794766655203362015666896734585),
                    std::complex<double>(18.120544172173429490704087812923, 5.4204235869878144503289374616496),
                    std::complex<double>(18.490141199340384891975196574175, 3.2548348766449188728805588482706),
                    std::complex<double>(18.673903714987714948137831741339, 1.0853727152211989477412078262598)
                };

                _residue_inf = 6.52606242395073216609341407232e+06;
            }
                  break;

            default: {
                std::cerr << "Error: the minor order of 32th-order IPF-PRAM method must be 8 !! \n";
                // RuntimeInformation::LogError("Position: void IpfPramSolver::CalcResiduesPoles; \n Error: the minor order of 32th-order IPF-PRAM method must be either 4 or 8 !!");
            }
                   break;
            }
        }
               break;

        case 48: {
            switch (MINOR_ORDER) {
            case 16: {
                _residues = {
                    std::complex<double>(66.970354207806730917746472111619, 33.085882130535162134156593476282),
                    std::complex<double>(81.944884359167433886220815117012, 17.681841799511539750929870266431),
                    std::complex<double>(92.27284637600642227891118687086, 2.1428314758051043362955739877301),
                    std::complex<double>(99.614424840437363374505995069078, -13.974856335115527491809200167586),
                    std::complex<double>(104.52058448308606031265599874563, -30.883926747214948265850169987144),
                    std::complex<double>(107.11900917515077472283145824826, -48.803743723623548729428550150017),
                    std::complex<double>(107.17284353427555676882654377297, -68.079341402195360598090051764673),
                    std::complex<double>(103.6517597589504140218428100027, -89.487576553120188134086191431648),
                    std::complex<double>(0, -27.848819652785689914669972405395),
                    std::complex<double>(0, -29.722431938406876838740718981203),
                    std::complex<double>(0, -31.876972733116498925400548731777),
                    std::complex<double>(0, -34.380029441533569787197533789292),
                    std::complex<double>(0, -37.322588844001770451220827102664),
                    std::complex<double>(0, -40.830228285502890427728931400369),
                    std::complex<double>(0, -45.08135431669163293659330516795),
                    std::complex<double>(0, -50.338301963183071695324322203336),
                    std::complex<double>(0, -57.003303913986738241982264812261),
                    std::complex<double>(0, -65.726116715499225118186938330207),
                    std::complex<double>(0, -77.629052184304011951186747847635),
                    std::complex<double>(0, -94.832046606338661533298951613979),
                    std::complex<double>(0, -121.87784511037333263178949986252),
                    std::complex<double>(0, -170.57762533959922433704347623647),
                    std::complex<double>(0, -284.23911224934968011040568355809),
                    std::complex<double>(0, -852.63206323866008046416497754505)
                };

                _poles = {
                    std::complex<double>(-20.415334484292075895957503964789, 53.99593932414685863539478425796),
                    std::complex<double>(-12.532755409033220327808957099614, 51.589604192952386104372247380539),
                    std::complex<double>(-6.5624782712902153476842129152373, 49.316998848474399616741278709945),
                    std::complex<double>(-1.6395386456417791172285254398134, 47.083807399104612872426519641001),
                    std::complex<double>(2.5748664677208408404487144375616, 44.86063510209574469965435948314),
                    std::complex<double>(6.2585829624309477970640766409355, 42.635516295139633453945468060203),
                    std::complex<double>(9.5197444778034849119359861739072, 40.403083877020324579650760895251),
                    std::complex<double>(12.43026626974239528100106213611, 38.160939786925954909460255955448),
                    std::complex<double>(15.040954926391404739518547303459, 35.908164599714767014201452765754),
                    std::complex<double>(17.389289058868392374006198615364, 33.644622420947161136997937926573),
                    std::complex<double>(19.503816734595369067921905404658, 31.37060750317470737526943429196),
                    std::complex<double>(21.406818989730454826349618880229, 29.086653392795744713187825751413),
                    std::complex<double>(23.116011785095083175007279875526, 26.793425401965735181285029239013),
                    std::complex<double>(24.645680662956471316615406784191, 24.491658312747133333234846393996),
                    std::complex<double>(26.007463494512089102229625634864, 22.182119751219279571578695466572),
                    std::complex<double>(27.210905629729082413128999936559, 19.865588647217181718898147689372),
                    std::complex<double>(28.263862513887334877989970746334, 17.542842806250618916296376241762),
                    std::complex<double>(29.172796836798126232847001165388, 15.214652104407450591689867501694),
                    std::complex<double>(29.943000667101314309363923673215, 12.881775209954092286894342046077),
                    std::complex<double>(30.578762790890893246183138128006, 10.544958542876781875906860151942),
                    std::complex<double>(31.083494955582070604775327465682, 8.2049366650222097071763783534089),
                    std::complex<double>(31.459826434882500249727929541004, 5.8624335871080518391006686597104),
                    std::complex<double>(31.709673423446690203680750531166, 3.5181646610363276366936365222801),
                    std::complex<double>(31.834287728092345118883736345296, 1.172838840005117410266661253415)
                };

                _residue_inf = 0.59332028180696474663102635488591;
            }
                   break;

            default: {
                std::cerr << "Error: the minor order of 48th-order IPF-PRAM method must be 16 !! \n";
                // RuntimeInformation::LogError("Position: void IpfPramSolver::CalcResiduesPoles; \n Error: the minor order of 48th-order IPF-PRAM method must be either 16 !!");
            }
                   break;
            }
        }
               break;

        case 64: {
            switch (MINOR_ORDER) {
            case 32: {
                _residues = {
                    std::complex<double>(111.67951638979344261570814368862, 51.385455632333670838443266864142),
                    std::complex<double>(128.93499687199647936504390379944, 33.214602045480550097278654981552),
                    std::complex<double>(141.76317666482498414348424487629, 16.039403264322539441644291961206),
                    std::complex<double>(151.98586275577407032168999922771, -0.97649149883691825961089446975354),
                    std::complex<double>(160.30489425780848151691882923548, -18.139882556089747618171639294761),
                    std::complex<double>(167.07694440123306487441889458631, -35.617120022390431668010926147006),
                    std::complex<double>(172.50546542105885839916537504354, -53.525568790761139364492029670321),
                    std::complex<double>(176.70960209244629526867627744385, -71.966367866349748090215254544038),
                    std::complex<double>(179.75327594895761868429487696996, -91.040107323912859146904947022284),
                    std::complex<double>(181.65720368508685575320519537507, -110.8568816271130836539074689618),
                    std::complex<double>(182.40091482142318196191425657571, -131.54514186773540773436973352972),
                    std::complex<double>(181.91552297000161710685818075754, -153.26209304867920669129033959558),
                    std::complex<double>(180.06235281320518893640243428926, -176.20901388236429028567214242135),
                    std::complex<double>(176.57996246342419232195132721919, -200.65885229098276561215692540102),
                    std::complex<double>(170.93486089295140415947440217126, -227.01987820163633200181279621099),
                    std::complex<double>(161.7182376249155024291422196087, -256.05994981495746837412785332499),
                    std::complex<double>(0, -24.133330174678789444636264103796),
                    std::complex<double>(0, -25.802481217476505894177801008291),
                    std::complex<double>(0, -27.717970177269877724663815792336),
                    std::complex<double>(0, -29.939092768600928566656394012029),
                    std::complex<double>(0, -32.545737099356359495989030114359),
                    std::complex<double>(0, -35.64819607035613991619896889097),
                    std::complex<double>(0, -39.403175480855305109679472373941),
                    std::complex<double>(0, -44.041100624056145955418740925488),
                    std::complex<double>(0, -49.915266988059246749190130138709),
                    std::complex<double>(0, -57.596359785292555050167010651435),
                    std::complex<double>(0, -68.070084718986097908542574125568),
                    std::complex<double>(0, -83.198305483597973381554790968655),
                    std::complex<double>(0, -106.97070169688692061151911927318),
                    std::complex<double>(0, -149.76040471243075038356376196437),
                    std::complex<double>(0, -249.6021750512315944233666829623),
                    std::complex<double>(0, -748.80870476559193358526203736532)
                };

                _poles = {
                    std::complex<double>(-19.326378328417174311105639151153, 87.204015730263863586017858396967),
                    std::complex<double>(-10.500411584685948998756297148744, 83.546875501994361392188813935213),
                    std::complex<double>(-3.6878857064120685614024740733091, 80.30420560265828030431205818352),
                    std::complex<double>(2.0261736729845799325149573354897, 77.252719050178981737335263187724),
                    std::complex<double>(6.9990684656214571519606574197488, 74.31183702019699746402691305267),
                    std::complex<double>(11.418202293988921239680781055455, 71.442276032832320180482325749955),
                    std::complex<double>(15.397632026267158389336104527071, 68.621748649449767664062215053),
                    std::complex<double>(19.013408149732517322095626021422, 65.836379055395374370405720621173),
                    std::complex<double>(22.319479936955594301606453470409, 63.076967733088448005648759389025),
                    std::complex<double>(25.355870695298764911623816074431, 60.337133159338428198710540969737),
                    std::complex<double>(28.153292315976092813682637496057, 57.61229582307771440177135671916),
                    std::complex<double>(30.735938571465171568176161781202, 54.899082241718600151070068584117),
                    std::complex<double>(33.123270331991049703951270261831, 52.194955704646370340309080860836),
                    std::complex<double>(35.331207165369130399385950933651, 49.497977262665210118028481248032),
                    std::complex<double>(37.372951471022606817914709223763, 46.806645382396531306121313507251),
                    std::complex<double>(39.259575560292766248509565206646, 44.119785112958077393725305499046),
                    std::complex<double>(41.000450410915662672079397917549, 41.436469511746935330866774729757),
                    std::complex<double>(42.603565492505142226141906532277, 38.755962714263354998186817369257),
                    std::complex<double>(44.07577169504004852692132715908, 36.077677896487890372194832894898),
                    std::complex<double>(45.422968731394944909025910572767, 33.401145710359164928173300606394),
                    std::complex<double>(46.650251629654385506945455738446, 30.725990225607042460717408397991),
                    std::complex<double>(47.762026529871881526376623982996, 28.051910341448298564684583607335),
                    std::complex<double>(48.762103060944555270116229999448, 25.378665241989615856759407003847),
                    std::complex<double>(49.653768566596028563720176366201, 22.706062878314617717031989750165),
                    std::complex<double>(50.439848051357433831625311693899, 20.0339507397450255534479102547),
                    std::complex<double>(51.122752724926477569671566307793, 17.362208370942111585030424301975),
                    std::complex<double>(51.704519306087918869822673443215, 14.690741228372236015288250750023),
                    std::complex<double>(52.186841720059609544088118873176, 12.019475567289573682416349413841),
                    std::complex<double>(52.571096428510864881471046517725, 9.3483541206788416208856710691446),
                    std::complex<double>(52.858362330558649129396644657006, 6.6773323824825090628596092592913),
                    std::complex<double>(53.049435938599895157626885669008, 4.0063753442643158187819256906766),
                    std::complex<double>(53.144842345525882885796444135442, 1.3354545608721807279187424262131)
                };

                _residue_inf = 4.82219923991115040585667182304943120892167097684066304e+05;
            }
                   break;

            default: {
                std::cerr << "Error: the minor order of 64th-order IPF-PRAM method must be 32 !! \n";
                // RuntimeInformation::LogError("Position: void IpfPramSolver::CalcResiduesPoles; \n Error: the minor order of 64th-order IPF-PRAM method must be 32 !!");
            }
                   break;
            }
        }
               break;

        default: {
            std::cerr << "Error: the order of IPF-PRAM method must be chosen from 16, 32, 48 and 64 !! \n";
            // RuntimeInformation::LogError("Position: void IpfPramSolver::CalcResiduesPoles; \n Error: the order of the IPF-PRAM method must be chosen from 16, 32, 48 and 64 !!");
        }
               break;
        }
    }
};

#endif