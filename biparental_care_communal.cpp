
///////////////////////////////
//          
//          Sex-biased dispersal and the evolution of sex-differences in care
//          when parental care is communal
//           
//          Individual-based simulations
//
//          This source code is used in the paper:
//          Kuijper, B. & Johnstone RA. (2015) 
//          How sex-biased dispersal affects conflict over parental investment/
//           
//           
//          (c) Bram Kuijper & Rufus A Johnstone
//          
//
//
//     This program is free software: you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
// 
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
// 
//     You should have received a copy of the GNU General Public License
//     along with this program.  If not, see <http://www.gnu.org/licenses/>. 
//
//
///////////////////////////////



#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <cmath>
#include <cassert>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "auxiliary.h"


//#define NDEBUG

using namespace std;

// random number generator 
// see http://www.gnu.org/software/gsl/manual/html_node/Random-Number-Generation.html#Random-Number-Generation 
gsl_rng_type const * T; // gnu scientific library rng type
gsl_rng *r; // gnu scientific rng 

const int Npatches = 200; 
const int Nfp = 4; // females per patch
const int Nmp = 2; // males per patch
const int numgen = 500000;
const int Clutch = 200;

// mutation rates
double mu = 0;
double sdmu = 0;

double df = 0; // female-specific dispersal
double dm = 0; // male-specific dispersal
double kf = 0; // cost of female care
double km = 0; // cost of male care
double init_um = 0; // initial level of male care
double init_uf = 0; // initial level of female care

// stats to count dispersers
int NdispF = 0;
int NdispM = 0;

// runtime for simulation
time_t total_time; 

// current generation
int generation = 0;

int seed = -1;

// skip the number of generations in output
// to prevent output files from becoming too large
int skip = 10;

// haploid individual
struct Individual
{
    // diploid locus of female parental care
    double uf[2];

    // diploid locus of male parental care
    double um[2];

    // the phenotypes of uf and um
    double uf_phen;
    double um_phen;
};

struct Patch
{
    Individual localsF[Nfp]; // all the local female breeders
    Individual localsM[Nmp]; // all the local male breeders

    // philopatric female offspring
    Individual philsF[Nfp * Clutch];     

    // philopatric male offspring
    Individual philsM[Nfp * Clutch];     

    // (note that dispersing offspring
    // ends up in global pool, see below)

    // total number of kids 
    int NkidsF; 
    int NkidsM; 

    // variable that allows for correct
    // rounding of the number of immigrants
    // per patch (see below)
    int immigrant_bonusF;
    int immigrant_bonusM;

};

// female and male mortality functions
double mortality_f(double const u)
{
    return(kf + (1.0-kf) * u * u); 
}

double mortality_m(double const u)
{
    return(km + (1.0-km) * u * u); 
}

// generate the complete population
Patch MetaPop[Npatches];
Individual DispersersF[Npatches * Nfp * Clutch];
Individual DispersersM[Npatches * Nfp * Clutch];

// give the outputfile a unique name
// by using the create_filename function (see 
// "bramauxiliary.h")
string filename("sim_biparental_care");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  // output file 

// initialize the command line arguments to vary 
// the parameters
void init_arguments(int argc, char *argv[])
{
    mu = atof(argv[1]); 
    sdmu = atof(argv[2]); 
    df = atof(argv[3]);
    dm = atof(argv[4]);
    kf = atof(argv[5]);
    km = atof(argv[6]);
    init_uf = atof(argv[7]);
    init_um = atof(argv[8]);
}

void init_pop()
{
    // start the time
    total_time = time(NULL);

    // obtain a seed from current nanosecond count
	seed = get_nanoseconds();
    
    // set up the random number generators
    // (from the gnu gsl library)
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed);

    // go through all patches
    for (int i = 0; i < Npatches; ++i)
    {
        // initialize genotypic and phenotypic values
        //
        // in females
        for (int j = 0; j < Nfp; ++j)
        {
            MetaPop[i].localsF[j].uf[0] = init_uf;
            MetaPop[i].localsF[j].uf[1] = init_uf;
            MetaPop[i].localsF[j].um[0] = init_um;
            MetaPop[i].localsF[j].um[1] = init_um;
            MetaPop[i].localsF[j].um_phen = init_um;
            MetaPop[i].localsF[j].uf_phen = init_uf;
        }
        // and in males
        for (int j = 0; j < Nmp; ++j)
        {
            MetaPop[i].localsM[j].uf[0] = init_uf;
            MetaPop[i].localsM[j].uf[1] = init_uf;
            MetaPop[i].localsM[j].um[0] = init_um;
            MetaPop[i].localsM[j].um[1] = init_um;
            MetaPop[i].localsM[j].um_phen = init_um;
            MetaPop[i].localsM[j].uf_phen = init_uf;
        }
    }
}

// mutate an allele
double Mut(double h)
{
    h+=gsl_rng_uniform(r) < mu ? gsl_ran_gaussian(r, sdmu) : 0;
    h = h < 0 ? 0 : h > 1.0 ? 1.0 : h;
    return(h);
}

// allocate a kid and give it genes
void Create_Kid(Individual &mother, Individual &father, Individual &Kid)
{
    // inherit (and potentially mutate) randomly sampled maternal allele for uf
    Kid.uf[0] = Mut(mother.uf[gsl_rng_uniform_int(r,2)]); 
    // inherit randomly sampled paternal allele for uf
    Kid.uf[1] = Mut(father.uf[gsl_rng_uniform_int(r,2)]);
    
    // inherit randomly sampled maternal allele for um
    Kid.um[0] = Mut(mother.um[gsl_rng_uniform_int(r,2)]);
    // inherit randomly sampled paternal allele for um
    Kid.um[1] = Mut(father.um[gsl_rng_uniform_int(r,2)]);

    // express phenotypes
    Kid.um_phen = .5 * (Kid.um[0] + Kid.um[1]);
    Kid.uf_phen = .5 * (Kid.uf[0] + Kid.uf[1]);
}

// mate and create kids across all patches...
void make_juveniles()
{
    // reset counters of dispersing
    // female and male juveniles
    NdispF = 0;
    NdispM = 0;

    // auxiliary variable storing randomly chosen father
    int random_father;
    // auxiliary variable storing randomly chosen mother 
    int random_mother;

    double uf, um;
    
    double skew_ratio = (double) Nmp / Nfp;

    // looping through the patches and creating offspring
    for (int i = 0; i < Npatches; ++i)
    {
        // reset philopatric juvenile counters
        MetaPop[i].NkidsF = 0;
        MetaPop[i].NkidsM = 0;

        // variable that takes into account
        // any beyond-the-minimum number of immigrants
        // (see later at the replace_adults function)
        MetaPop[i].immigrant_bonusF = 0;
        MetaPop[i].immigrant_bonusM = 0;

        // let each local female produce offspring
        for (int j = 0; j < Nfp; ++j)
        {
            // create kids 
            for (int k = 0; k < Clutch; ++k)
            {
                // select randomly chosen local father
                // to provide care
                random_father = gsl_rng_uniform_int(r, Nmp);
                // select randomly chosen local mother
                // that is going to perform care for this
                // offspring
                random_mother = gsl_rng_uniform_int(r, Nfp);

                // get social mother's value of parental care
                // (not necessarily the same as the biological mother)
                uf = MetaPop[i].localsF[random_mother].uf_phen;

                // get father's value of parental care
                um = MetaPop[i].localsM[random_father].um_phen;

                // calculate offspring survival, b
                if (gsl_rng_uniform(r) > 1.0 - exp(-(uf + skew_ratio * um)))
                {
                    // death, continue and make next offspring
                    continue;
                }

                // create an offspring
                Individual Kid; 

                // father siring the offspring might be different 
                // from father giving care
                random_father = gsl_rng_uniform_int(r, Nmp);

                Create_Kid(MetaPop[i].localsF[j],MetaPop[i].localsM[random_father],Kid);

                // female or male
                if (gsl_rng_uniform(r) < 0.5)
                {
                    // disperse or not
                    if (gsl_rng_uniform(r) < df)
                    {
                        DispersersF[NdispF++] = Kid;
                    }
                    else
                    {
                        MetaPop[i].philsF[MetaPop[i].NkidsF++] = Kid;
                    }
                }
                else
                {
                    if (gsl_rng_uniform(r) < dm)
                    {
                        DispersersM[NdispM++] = Kid;
                    }
                    else
                    {
                        MetaPop[i].philsM[MetaPop[i].NkidsM++] = Kid;
                    }
                }
            } //clutch
        }//Npp
    }//Npatches

    assert(NdispF < Npatches * Nfp * Clutch);
    assert(NdispM < Npatches * Nfp * Clutch);
}

// replacement of adults with juveniles
void replace_adults()
{
    // we first need to distribute Ndisp individuals over the various 
    // patches. As it is unlikely that Ndisp/Npatches is exactly an integer
    // we need to some rounding
    
    // first, I assign to each patch just the minimum number of dispersers
    // This is simply floor(Ndisp/Npatches)
    int dispersers_per_patchF = floor((double) NdispF / Npatches);
    int dispersers_per_patchM = floor((double) NdispM / Npatches);

    // Then a certain randomly selected number of patches may 
    // receive dispersers beyond that minimum until the pool of dispersers
    // is depleted
    for (int i = 0; i < NdispF - Npatches * dispersers_per_patchF; ++i)
    {
        // randomly picked patch receives additional immigrant
        MetaPop[gsl_rng_uniform_int(r,Npatches)].immigrant_bonusF++;
    }
    for (int i = 0; i < NdispM - Npatches * dispersers_per_patchM; ++i)
    {
        // randomly picked patch receives additional immigrant
        MetaPop[gsl_rng_uniform_int(r,Npatches)].immigrant_bonusM++;
    }
    
    // now replace local breeders on each patch
    for (int i = 0; i < Npatches; ++i)
    {
        if (dispersers_per_patchF + MetaPop[i].immigrant_bonusF > 0)
        {
            assert(NdispF > 0);
        }
        if (dispersers_per_patchM + MetaPop[i].immigrant_bonusM > 0)
        {
            assert(NdispM > 0);
        }

        int arriving_immigrantsF = dispersers_per_patchF + MetaPop[i].immigrant_bonusF;
        int arriving_immigrantsM = dispersers_per_patchM + MetaPop[i].immigrant_bonusM;

        for (int j = 0; j < Nfp; ++j)
        {
            assert((double) arriving_immigrantsF / (arriving_immigrantsF + MetaPop[i].NkidsF) >= 0.0 
                    && (double) arriving_immigrantsF / (arriving_immigrantsF + MetaPop[i].NkidsF) <= 1.0);
            
            // female dies with a certain probability
            if (gsl_rng_uniform(r) < mortality_f(0.5 * (MetaPop[i].localsF[j].uf[0] + MetaPop[i].localsF[j].uf[1])))
            {
                // replace her either with an immigrant or philopatrically born juvenile
                if (gsl_rng_uniform(r) < (double) arriving_immigrantsF / (arriving_immigrantsF + MetaPop[i].NkidsF))
                {
                    int rand_disp = gsl_rng_uniform_int(r,NdispF);
                    MetaPop[i].localsF[j] = DispersersF[rand_disp];
                    DispersersF[rand_disp] = DispersersF[NdispF-1];
                    --NdispF;
                    --arriving_immigrantsF;
                }
                else
                {
                    int rand_phil = gsl_rng_uniform_int(r,MetaPop[i].NkidsF);
                    MetaPop[i].localsF[j] = MetaPop[i].philsF[rand_phil];
                    MetaPop[i].philsF[rand_phil] = MetaPop[i].philsF[MetaPop[i].NkidsF-1];
                    --MetaPop[i].NkidsF;
                }
            }
        }
        
        for (int j = 0; j < Nmp; ++j)
        {
            assert((double) arriving_immigrantsM / (arriving_immigrantsM + MetaPop[i].NkidsM) >= 0.0 
                    && (double) arriving_immigrantsM / (arriving_immigrantsM + MetaPop[i].NkidsM) <= 1.0);

            // male dies 
            if (gsl_rng_uniform(r) < mortality_m(0.5 * (MetaPop[i].localsM[j].um[0] + MetaPop[i].localsM[j].um[1])))
            {
                if (gsl_rng_uniform(r) < (double) arriving_immigrantsM / (arriving_immigrantsM + MetaPop[i].NkidsM))
                {
                    int rand_disp = gsl_rng_uniform_int(r,NdispM);
                    MetaPop[i].localsM[j] = DispersersM[rand_disp];
                    DispersersM[rand_disp] = DispersersM[NdispM-1];
                    --NdispM;
                    --arriving_immigrantsM;
                }
                else
                {
                    int rand_phil = gsl_rng_uniform_int(r,MetaPop[i].NkidsM);
                    MetaPop[i].localsM[j] = MetaPop[i].philsM[rand_phil];
                    MetaPop[i].philsM[rand_phil] = MetaPop[i].philsM[MetaPop[i].NkidsM-1];
                    --MetaPop[i].NkidsM;
                }
            }
        }

        // remove all remaining female immigrants from the global dispersal pool
        for (int j = 0; j < arriving_immigrantsF; ++j)
        {
            int rand_disp = gsl_rng_uniform_int(r,NdispF);
            DispersersF[rand_disp] = DispersersF[NdispF-1];
            --NdispF;
        }
        // remove all remaining male immigrants from the global dispersal pool
        for (int j = 0; j < arriving_immigrantsM; ++j)
        {
            int rand_disp = gsl_rng_uniform_int(r,NdispM);
            DispersersM[rand_disp] = DispersersM[NdispM-1];
            --NdispM;
        }
    }
    assert(NdispF==0);
    assert(NdispM==0);
}

void write_data_headers()
{
    DataFile << "generation;uf;um;clutch;varuf;varum;" << endl;
}

void write_data()
{
    double meanuf = 0;
    double meanum = 0;
    double ssuf = 0;
    double ssum = 0;

    double uf = 0;
    double um = 0;

    for (int i = 0; i < Npatches; ++i)
    {
        for (int j = 0; j < Nfp; ++j)
        {
            uf = 0.5 * (MetaPop[i].localsF[j].uf[0] + MetaPop[i].localsF[j].uf[1]);
            um = 0.5 * (MetaPop[i].localsF[j].um[0] + MetaPop[i].localsF[j].um[1]);
            meanuf += uf;
            meanum += um;
            ssuf += uf * uf;
            ssum += um * um;
        }
        
        for (int j = 0; j < Nmp; ++j)
        {
            uf = 0.5 * (MetaPop[i].localsM[j].uf[0] + MetaPop[i].localsM[j].uf[1]);
            um = 0.5 * (MetaPop[i].localsM[j].um[0] + MetaPop[i].localsM[j].um[1]);
            meanuf += uf;
            meanum += um;
            ssuf += uf * uf;
            ssum += um * um;
        }
    }
    
    meanuf /= Npatches * (Nfp + Nmp);
    meanum /= Npatches * (Nfp + Nmp);
        
    DataFile << generation << ";" << meanuf << ";" << meanum << ";" 
        << (2 * (meanuf + meanum) - (meanuf + meanum)*(meanuf+meanum))*Clutch << ";" <<
                                ssuf / (Npatches * (Nfp + Nmp))  - meanuf * meanuf << ";" 
                                << ssum / (Npatches * (Nfp + Nmp))  - meanum * meanum << ";" << endl;

}

void write_parameters()
{
    DataFile << endl << endl << "patch;" << Npatches << endl
                << "nfp;" << Nfp << endl
                << "nmp;" << Nmp << endl
                << "numgen;" << numgen << endl
                << "mu;" << mu << endl
                << "sdmu;" << sdmu << endl
                << "kf;" << kf << endl
                << "km;" << km << endl
                << "df;" << df << endl
                << "dm;" << dm << endl
                << "runtime;" << total_time << endl;
}

// the central part of the code
int main(int argc, char * argv[])
{
    init_arguments(argc,argv);
    init_pop();

    write_data_headers();

    for (generation = 0; generation < numgen; ++generation)
    {
        make_juveniles();

        if (generation % skip == 0)
        {
            write_data();
        }

        replace_adults();
    }

    write_data();
    write_parameters();
}
