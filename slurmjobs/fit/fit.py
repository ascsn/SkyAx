# this is a python code that calculates the half-lives
# by calling the F77 solver
# the optimizer is called from the PETSc/TAO toolkit

import numpy as np
import os
import sys
import subprocess
import time
from petsc4py import PETSc



def calculate_response_K(nucleus,V0_pp,g_0,Jpi,k_val,types,recalc):
    '''
    Describe the module here

    Jpi can have the following values:
    1+ : GT response
    2- : FF 2m response
    1- : FF 1m response
    0- : FF 0m response

    This module calculates for only one K
    given by k_val

    types:
    This module is now modified, it has to handle both b- and EC/b+.
    Types can have values HL- and HL+, b+ always goes together with EC,
    so two calculations are needed. Therefore, I change the structure of 
    subdirectories in the following way:

    Fe50/b+/GT_1p/K0/... or Fe50/b-/GT_1p/K0 ...

    All the scripts have to be modified according to this.

    If recalc == True then the ..._fix.sh scripts will be run
    Note that these are used to avoid problems with scavenger queue

    NOTE: assumes that qrpa.wel file is already prepared

    '''

    # prepare parent directory
    folder_path='/mnt/gs21/scratch/ravlican/fit_me_beta_pm/chi2_fit/'
    os.makedirs(folder_path+nucleus+'/'+types, exist_ok=True)

    # NOTE: assumes that the ground state is precalculated, so it can copy qrpa.wel files
    # now move the qrpa.wel file into the directory for that nucleus
    subprocess.run(['cp /mnt/gs21/scratch/ravlican/dirhbz_blocking_QRPA/fit/'+nucleus+'/qrpa.wel /mnt/gs21/scratch/ravlican/fit_me_beta_pm/chi2_fit/'+nucleus+'/'+types], shell=True)

    # now go to the subfolder with the nucleus
    os.chdir('/mnt/gs21/scratch/ravlican/fit_me_beta_pm/chi2_fit/'+nucleus+'/'+types)

    # create a subfolder for the GT mode
    if (Jpi == '1+'):
        folder_path='/mnt/gs21/scratch/ravlican/fit_me_beta_pm/chi2_fit/'+nucleus+'/'+types+'/GT_1p'
    elif (Jpi == '1-'):
        folder_path='/mnt/gs21/scratch/ravlican/fit_me_beta_pm/chi2_fit/'+nucleus+'/'+types+'/FF_1m'
    elif (Jpi == '0-'):
        folder_path='/mnt/gs21/scratch/ravlican/fit_me_beta_pm/chi2_fit/'+nucleus+'/'+types+'/FF_0m'
    elif (Jpi == '2-'):
        folder_path='/mnt/gs21/scratch/ravlican/fit_me_beta_pm/chi2_fit/'+nucleus+'/'+types+'/FF_2m'
    else:
        print('Transition not implemented !', Jpi)
        sys.exit(-1)

    os.makedirs(folder_path, exist_ok=True)
    os.chdir(folder_path)


    os.makedirs('K'+str(k_val), exist_ok=True)

    # now copy all the files you need to run the RPA code
    subprocess.run(['cp ../../../../run K'+str(k_val)+'/.'], shell=True)
    subprocess.run(['cp ../../../../rpa.dat K'+str(k_val)+'/.'], shell=True)
    subprocess.run(['cp ../../../../rpa.par K'+str(k_val)+'/.'], shell=True)
    subprocess.run(['cp ../../../../*.f K'+str(k_val)+'/.'], shell=True)
    subprocess.run(['cp ../../../../makefile K'+str(k_val)+'/.'], shell=True)


    if (Jpi == '1+'):
        subprocess.run(['cp ../../../../script_GT_1p_run.sh K'+str(k_val)+'/.'], shell=True)
        subprocess.run(['cp ../../../../script_GT_1p_run_fix.sh K'+str(k_val)+'/.'], shell=True)
    elif (Jpi == '1-'):
        subprocess.run(['cp ../../../../script_FF_1m_run.sh K'+str(k_val)+'/.'], shell=True)
        subprocess.run(['cp ../../../../script_FF_1m_run_fix.sh K'+str(k_val)+'/.'], shell=True)
    elif (Jpi == '0-'):
        subprocess.run(['cp ../../../../script_FF_0m_run.sh K'+str(k_val)+'/.'], shell=True)
        subprocess.run(['cp ../../../../script_FF_0m_run_fix.sh K'+str(k_val)+'/.'], shell=True)
    else:
        subprocess.run(['cp ../../../../script_FF_2m_run.sh K'+str(k_val)+'/.'], shell=True)
        subprocess.run(['cp ../../../../script_FF_2m_run_fix.sh K'+str(k_val)+'/.'], shell=True)

    # modify the rpa.dat file
    file_path= 'K'+str(k_val)+'/rpa.dat'
    with open(file_path, 'r') as file:
        lines = file.readlines()

        # Modify the specific line
    if (Jpi == '1+'):
        new_line_content = 'J  Pi    =    1   +1              ! ang. mom. and parity           |\n'
    elif (Jpi == '1-'):
        new_line_content = 'J  Pi    =    1   -1              ! ang. mom. and parity           |\n'
    elif (Jpi == '0-'):
        new_line_content = 'J  Pi    =    0   -1              ! ang. mom. and parity           |\n'
    else:
        new_line_content = 'J  Pi    =    2   -1              ! ang. mom. and parity           |\n'

    lines[2 - 1] = new_line_content


    new_line_content = 'K        =    {}                   ! ang. mom. projection           |\n'.format(k_val)
    lines[3 - 1] = new_line_content

    if (types == 'b-'):
        new_line_content = 'iso      =    0                   ! 0 (b-) 1 (b+)|\n'
    elif (types == 'b+'):
        new_line_content = 'iso      =    1                   ! 0 (b-) 1 (b+)|\n'
    elif (types == 'EC'):
        new_line_content = 'iso      =    1                   ! 0 (b-) 1 (b+)|\n'
    else:
        print('Wrong type value !')
        sys.exit(-1)

    lines[4 - 1] = new_line_content


    if (types == 'b-'):
        new_line_content = 'beta_dec =    1                   | 1b-, 2b+, 3EC                  |\n'
    elif (types == 'b+'):
        new_line_content = 'beta_dec =    2                   | 1b-, 2b+, 3EC                  |\n'
    elif (types == 'EC'):
        new_line_content = 'beta_dec =    3                   | 1b-, 2b+, 3EC                  |\n'
    else:
        print('Wrong type value !')
        sys.exit(-1)

    lines[9 - 1] = new_line_content

    new_line_content = 'V0_pp    =  {}                                                   |\n'.format(round(V0_pp,5))
    lines[15 - 1] = new_line_content

    new_line_content = 'deltafac =  {}                                                   |\n'.format(round(g_0,5))
    lines[17 - 1] = new_line_content

    if (Jpi == '1+'):
        new_line_content = 'N_oper   =  1                     ! number of operators to consdr  |\n'
    elif (Jpi == '1-'):
        new_line_content = 'N_oper   =  5                     ! number of operators to consdr  |\n'
    elif (Jpi == '0-'):
        new_line_content = 'N_oper   =  3                     ! number of operators to consdr  |\n'
    else:
        new_line_content = 'N_oper   =  1                     ! number of operators to consdr  |\n'

    lines[19 - 1] = new_line_content

    if (Jpi == '1+'):
        new_line_content = 'GAMT\n'
        lines[20 - 1] = new_line_content
    elif (Jpi == '1-'):
        new_line_content = 'NSDR\n'
        lines[20 - 1] = new_line_content
        new_line_content = 'RSDR\n'
        lines[21 - 1] = new_line_content
        new_line_content = 'NDIP\n'
        lines[22 - 1] = new_line_content
        new_line_content = 'RDIP\n'
        lines[23 - 1] = new_line_content
        new_line_content = 'ALPH\n'
        lines[24 - 1] = new_line_content
    elif (Jpi == '0-'):
        new_line_content = 'NSDR\n'
        lines[20 - 1] = new_line_content
        new_line_content = 'RSDR\n'
        lines[21 - 1] = new_line_content
        new_line_content = 'GAM5\n'
        lines[22 - 1] = new_line_content
    else:
        new_line_content = 'NSDR\n'
        lines[20 - 1] = new_line_content




    # Write the updated content back to the file
    with open(file_path, 'w') as file:
        file.writelines(lines)

    # now modify the script for the specific problem
    # Read the content of the file
    if recalc:
        if (Jpi == '1+'):
            file_path= 'K'+str(k_val)+'/script_GT_1p_run_fix.sh'
        elif (Jpi == '1-'):
            file_path= 'K'+str(k_val)+'/script_FF_1m_run_fix.sh'
        elif (Jpi == '0-'):
            file_path= 'K'+str(k_val)+'/script_FF_0m_run_fix.sh'
        else:
            file_path= 'K'+str(k_val)+'/script_FF_2m_run_fix.sh'
    else:
        if (Jpi == '1+'):
            file_path= 'K'+str(k_val)+'/script_GT_1p_run.sh'
        elif (Jpi == '1-'):
            file_path= 'K'+str(k_val)+'/script_FF_1m_run.sh'
        elif (Jpi == '0-'):
            file_path= 'K'+str(k_val)+'/script_FF_0m_run.sh'
        else:
            file_path= 'K'+str(k_val)+'/script_FF_2m_run.sh'

    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Modify the specific line
    new_line_content = '#SBATCH --job-name='+nucleus+'_K'+str(k_val)+'_'+types+'\n'
    lines[2 - 1] = new_line_content

    result = subprocess.run(['pwd'], stdout=subprocess.PIPE)
    current_directory = result.stdout.strip().decode('utf-8')
    print('Current directory, ', current_directory, flush=True)
    new_line_content = 'dir='+current_directory+'/K'+str(k_val)+'\n'
    
    lines[12 - 1] = new_line_content  

    new_line_content = 'K_val={}\n'.format(k_val)
    lines[15 - 1] = new_line_content 

    new_line_content = 'nucleus={}\n'.format(nucleus)
    lines[16 - 1] = new_line_content 

    

    # Write the updated content back to the file
    with open(file_path, 'w') as file:
        file.writelines(lines)

    # at this point you can sbatch the files
    if recalc:
        if (Jpi == '1+'):
            subprocess.run(['sbatch K'+str(k_val)+'/script_GT_1p_run_fix.sh'], shell=True)
        elif (Jpi == '1-'):
            subprocess.run(['sbatch K'+str(k_val)+'/script_FF_1m_run_fix.sh'], shell=True)
        elif (Jpi == '0-'):
            subprocess.run(['sbatch K'+str(k_val)+'/script_FF_0m_run_fix.sh'], shell=True)
        else:
            subprocess.run(['sbatch K'+str(k_val)+'/script_FF_2m_run_fix.sh'], shell=True)
    else:
        if (Jpi == '1+'):
            subprocess.run(['sbatch K'+str(k_val)+'/script_GT_1p_run.sh'], shell=True)
        elif (Jpi == '1-'):
            subprocess.run(['sbatch K'+str(k_val)+'/script_FF_1m_run.sh'], shell=True)
        elif (Jpi == '0-'):
            subprocess.run(['sbatch K'+str(k_val)+'/script_FF_0m_run.sh'], shell=True)
        else:
            subprocess.run(['sbatch K'+str(k_val)+'/script_FF_2m_run.sh'], shell=True)

    
    print('Finished submitting response ...', flush=True)

    return




def block_jobs():

    # now check for the processes to finish before moving on with the script
    number_of_jobs=1
    while number_of_jobs:
        # Check the number of jobs in the queue using squeue
        squeue_command = 'squeue -u ravlican'  # Replace <username> with your actual username
        result = subprocess.run(squeue_command, shell=True, stdout=subprocess.PIPE)

        # Decode the byte output to string
        output = result.stdout.decode('utf-8') if result.stdout else ''

        # Process the output to count the number of jobs in the queue
        output_lines = output.split('\n')
        # Exclude the header line and count the remaining lines (each line represents a job)
        number_of_jobs = len(output_lines) - 2 if output_lines else 0 # was -2 !

        print(f"Number of jobs in the queue: {number_of_jobs}", flush=True)
        
        # Wait for a specified time interval before checking again
        time.sleep(10)  # Adjust the interval based on how frequently you want to check

    return



def check_jobs(param,user_data):
    '''
    Checks whether all the calculations
    have converged to something

    '''

    data, nucleus, types = user_data # unpack tuple

    safety_counter = 0
    total_recalc = 1
    while (total_recalc > 0):
        '''
        This will loop until the condition is meet that all files have
        a required number of lines !

        '''
        recalc_bm = 0 # how many processes have to be recalculated
        recalc_bp = 0
        recalc_EC = 0

        for nuc in range(len(nucleus)):

            for jpi in ['1+','1-','0-','2-']:

                # 1- mode
                if (jpi == '1-'):
                    ops=['NSDRNSDR', 'NSDRRSDR', 'NSDRNDIP', 'NSDRRDIP', 'NSDRALPH',
                    'RSDRRSDR', 'RSDRNDIP', 'RSDRRDIP', 'RSDRALPH', 'NDIPNDIP',
                    'NDIPRDIP','NDIPALPH','RDIPRDIP','RDIPALPH','ALPHALPH']
                    label = '1m'
                    label_folder = 'FF_1m'
                    kmax = 2
                elif (jpi == '1+'):
                    ops=['GAMTGAMT']
                    label = '1p'
                    label_folder = 'GT_1p'
                    kmax = 2
                elif (jpi == '0-'):
                    ops=['NSDRNSDR','NSDRRSDR','RSDRRSDR','GAM5GAM5','RSDRGAM5','NSDRGAM5']
                    label = '0m'
                    label_folder = 'FF_0m'
                    kmax = 1
                else:
                    ops=['NSDRNSDR']
                    label = '2m'
                    label_folder = 'FF_2m'
                    kmax = 3

                if (types[nuc] == 'HL-'):

                    folder = '/mnt/gs21/scratch/ravlican/fit_me_beta_pm/chi2_fit/'+nucleus[nuc]+'/b-/'+label_folder

                    for k_val in range(kmax):
                        problem = False
                        for o in range(len(ops)):
                            file = folder+'/K'+str(k_val)+'/beta_points_'+label+'_K'+str(k_val)+'_'+ops[o]+'_'+nucleus[nuc]+'.dat'
                            try:
                                test_data = np.loadtxt(file, skiprows=1)
                                if (len(test_data) != 60 ):
                                    print('Problem with file, ', file)
                                    problem = True
                            except FileNotFoundError:
                                print('File not created: ', file)
                                problem = True

                        # if there is a problem, restart the calculation for that K value
                        if (problem):
                            get_me_response_K(param,nucleus[nuc],jpi,k_val,'b-',True)
                            recalc_bm += 1

                if (types[nuc] == 'HL+'):

                    folder = '/mnt/gs21/scratch/ravlican/fit_me_beta_pm/chi2_fit/'+nucleus[nuc]+'/b+/'+label_folder

                    
                    
                    for k_val in range(kmax):
                        problem = False
                        for o in range(len(ops)):
                            file = folder+'/K'+str(k_val)+'/beta_points_'+label+'_K'+str(k_val)+'_'+ops[o]+'_'+nucleus[nuc]+'.dat'
                            try:
                                test_data = np.loadtxt(file, skiprows=1)
                                if (len(test_data) != 60 ):
                                    print('Problem with file, ', file)
                                    problem = True
                            except FileNotFoundError:
                                print('File not created: ', file)
                                problem = True

                        # if there is a problem, restart the calculation for that K value
                        if (problem):
                            get_me_response_K(param,nucleus[nuc],jpi,k_val,'b+',True)
                            recalc_bp += 1



                    folder = '/mnt/gs21/scratch/ravlican/fit_me_beta_pm/chi2_fit/'+nucleus[nuc]+'/EC/'+label_folder

                    
                    
                    for k_val in range(kmax):
                        problem = False
                        for o in range(len(ops)):
                            file = folder+'/K'+str(k_val)+'/beta_points_'+label+'_K'+str(k_val)+'_'+ops[o]+'_'+nucleus[nuc]+'.dat'
                            try:
                                test_data = np.loadtxt(file, skiprows=1)
                                if (len(test_data) != 60 ):
                                    print('Problem with file, ', file)
                                    problem = True
                            except FileNotFoundError:
                                print('File not created: ', file)
                                problem = True

                        # if there is a problem, restart the calculation for that K value
                        if (problem):
                            get_me_response_K(param,nucleus[nuc],jpi,k_val,'EC',True)
                            recalc_EC += 1


        # wait until all the corrected jobs are finished
        total_recalc = recalc_bm + recalc_bp + recalc_EC
        print('Finished checking for problems before HL calculation!')
        print('Number of calculations that have to be restarted: ', total_recalc)
        print('Iteration ?', safety_counter)
        safety_counter += 1
        # increased safety limit because of stuff
        if (safety_counter > 50):
            print('Safety counter reached > 50 !')
            sys.exit(-1)

        block_jobs()



    print('All files have been checked and I am moving on to calculate HLs ...')

    return




#  define the objective function
def residual_routine(tao,solution,res,*user_data): 
    '''

    definition of the objective function for the
    minimization

    tao - TAO object
    solution - Vec() of parameters:
    solution[0] - isoscalar pairing strength 
    solution[1] - axial coupling g_A
    solution[2] - Landau-Migdal coupling
    res - Vec() of residuals
    user_data is tuple of following:
    data - list of exp data
    weights - list of errors
    nucleus - list of nuclei to calculate
    deformations - list of deformations for cnstr calc

    Remember that EC calculations need input from EC_files

    '''
    global iteration_num 

    param = solution.getArray() # turns Vec() to array

    data, nucleus, types = user_data # unpack tuple

    n_d = len(nucleus) # number of nuclei
    n_x = len(param) # number of parameters

    # I have to calculate all half-lives for a given
    # value of the parameter
    x = [] # list of calculated total half-lives
    x_GT = [] # list of calculated GT half-lives
    x_FF_1m = [] 
    x_FF_0m = []
    x_FF_2m = []
    x_GT_EC = [] # list of calculated GT half-lives
    x_FF_1m_EC = [] 
    x_FF_0m_EC = []
    x_FF_2m_EC = []
    E_GTR = [] # list of resonance energies
    for nuc in range(len(nucleus)):

        if (types[nuc] == 'HL-'):
            #continue
            # gets to the queue
            for k in range(2):
                get_me_response_K(param,nucleus[nuc],'1+',k,'b-',False)
            for k in range(2):
                get_me_response_K(param,nucleus[nuc],'1-',k,'b-',False)
            
            get_me_response_K(param,nucleus[nuc],'0-',0,'b-',False)

            for k in range(3):
               get_me_response_K(param,nucleus[nuc],'2-',k,'b-',False)

        elif (types[nuc] == 'HL+'):
            #continue
            for k in range(2):
                get_me_response_K(param,nucleus[nuc],'1+',k,'b+',False)
            for k in range(2):
                get_me_response_K(param,nucleus[nuc],'1-',k,'b+',False)
            
            get_me_response_K(param,nucleus[nuc],'0-',0,'b+',False)

            for k in range(3):
               get_me_response_K(param,nucleus[nuc],'2-',k,'b+',False)

            for k in range(2):
                get_me_response_K(param,nucleus[nuc],'1+',k,'EC',False)
            for k in range(2):
                get_me_response_K(param,nucleus[nuc],'1-',k,'EC',False)
            
            get_me_response_K(param,nucleus[nuc],'0-',0,'EC',False)

            for k in range(3):
               get_me_response_K(param,nucleus[nuc],'2-',k,'EC',False)


        else:
            # calculates straight away
            E_GTR.append(calculate_GT_resonance(param,nucleus[nuc]))

    #block untill all processes are finished
    block_jobs()

    #check that all jobs have converged before moving on
    check_jobs(param,user_data)


    # now calculate the half_lives
    for nuc in range(len(nucleus)):

        if (types[nuc] == 'HL-'):
            fac1 = GT_rate_module.calculate_HL_1p(nucleus[nuc],param[1])
            fac2 = FF_1m_rate_module.calculate_HL_1m(nucleus[nuc],param[1])
            fac3 = FF_0m_rate_module.calculate_HL_0m(nucleus[nuc],param[1])
            fac4 = FF_2m_rate_module.calculate_HL_2m(nucleus[nuc],param[1])
            fac5 = 0.0
            fac6 = 0.0
            fac7 = 0.0
            fac8 = 0.0
        if (types[nuc] == 'HL+'):
            fac1 = GT_rate_module_Bp.calculate_HL_1p(nucleus[nuc],param[1])
            fac2 = FF_1m_rate_module_Bp.calculate_HL_1m(nucleus[nuc],param[1])
            fac3 = FF_0m_rate_module_Bp.calculate_HL_0m(nucleus[nuc],param[1])
            fac4 = FF_2m_rate_module_Bp.calculate_HL_2m(nucleus[nuc],param[1])
            fac5 = GT_rate_module_EC.calculate_HL_1p(nucleus[nuc],param[1])
            fac6 = FF_1m_rate_module_EC.calculate_HL_1m(nucleus[nuc],param[1])
            fac7 = FF_0m_rate_module_EC.calculate_HL_0m(nucleus[nuc],param[1])
            fac8 = FF_2m_rate_module_EC.calculate_HL_2m(nucleus[nuc],param[1])

        # now here correct for negative half-lives
        '''
            this can happen in the instances where V0_pp is too strong
            just assign 1e-10 to the rate in that case so that it is
            negligible

        '''

        if (fac1 <= 0.0):
            fac1 = 1e10
        if (fac2 <= 0.0):
            fac2 = 1e10
        if (fac3 <= 0.0):
            fac3 = 1e10
        if (fac4 <= 0.0):
            fac4 = 1e10
        if (fac5 <= 0.0):
            fac5 = 1e10
        if (fac6 <= 0.0):
            fac6 = 1e10
        if (fac7 <= 0.0):
            fac7 = 1e10
        if (fac8 <= 0.0):
            fac8 = 1e10


        # append the partial half-lives
        if (types[nuc] == 'HL-'):
            x_GT.append(fac1)
            x_FF_1m.append(fac2)
            x_FF_0m.append(fac3)
            x_FF_2m.append(fac4)
            x_GT_EC.append(0.0)
            x_FF_1m_EC.append(0.0)
            x_FF_0m_EC.append(0.0)
            x_FF_2m_EC.append(0.0)

            # append the total half-life
            x.append(1./(1./fac1 + 1./fac2 + 1./fac3 + 1./fac4))
        
        elif (types[nuc] == 'HL+'):
            x_GT.append(fac1)
            x_FF_1m.append(fac2)
            x_FF_0m.append(fac3)
            x_FF_2m.append(fac4)
            x_GT_EC.append(fac5)
            x_FF_1m_EC.append(fac6)
            x_FF_0m_EC.append(fac7)
            x_FF_2m_EC.append(fac8)

            # append the total half-life
            x.append(1./(1./fac1 + 1./fac2 + 1./fac3 + 1./fac4 + \
            1./fac5 + 1./fac6 + 1./fac7 + 1./fac8))

        


    # append the resonances at the end of the list
    for ee in range(len(E_GTR)):
        x.append(E_GTR[ee])
        x_GT.append(0.0)
        x_FF_1m.append(0.0)
        x_FF_0m.append(0.0)
        x_FF_2m.append(0.0)

    
    
    suma = 0.0
    '''
    The residual function is calculated as follows:

    for types[i]='HL' we calculate the differences between the log10
    and weight is set to 0.5

    for types[i]='GTR' we just calculate the energy difference
    and error is set to 0.3, but it can be manually changed ...

    change w_GTR to whatever you want

    '''
    for i in range(n_d):
        
        if (types[i] == 'HL-' or types[i] == 'HL+'):
            res[i] = ((np.log10(x[i])-np.log10(data[i]))/0.5)/np.sqrt(n_d - n_x)
            suma += ((np.log10(x[i])-np.log10(data[i]))/0.5)**2
        else:
            res[i] = ((x[i]-data[i])/0.30)/np.sqrt(n_d - n_x)
            suma += ((x[i]-data[i])/0.30)**2

    
    # increase the number of iterations by 1
    iteration_num += 1

    # write down the result for each cost function
    file_path = '/mnt/gs21/scratch/ravlican/fit_me_beta_pm/chi2_fit/chi2.txt'

    # Open the file in write mode ('w')
    with open(file_path, 'a') as file:
        if (iteration_num == 1):
            file.write("# starting iterations ... \n")
        file.write("{} {} {} {} {:e}\n".format(iteration_num,param[0],param[1],param[2],suma/(n_d - n_x) ))

    # before returning I want to make a progress file
    # it will save HL at each iterations for all considered nuclei
    file_path = '/mnt/gs21/scratch/ravlican/fit_me_beta_pm/chi2_fit/chi2_iteration_'+str(iteration_num)+'.txt'

    # Open the file in write mode ('w')
    with open(file_path, 'w') as file:
        # write the header to the file
        file.write("# {} {} {} {} {:e}\n".format(iteration_num,param[0],param[1],param[2],suma/(n_d - n_x) ))
        # write calculated HL for each nucleus
        for nuc in range(len(nucleus)):
            if (types[nuc] == 'HL-'):
                file.write("{} {} {} {} {:e} {:e} {:e} {:e} {:e} {:e} {:e} {:e} {:e} {:e}\n".format(param[0],param[1],param[2],\
                nucleus[nuc], x[nuc], data[nuc], x_GT[nuc], x_FF_1m[nuc], x_FF_0m[nuc], x_FF_2m[nuc],\
                 0.0, 0.0, 0.0, 0.0 ))
            elif (types[nuc] == 'HL+'):
                file.write("{} {} {} {} {:e} {:e} {:e} {:e} {:e} {:e} {:e} {:e} {:e} {:e}\n".format(param[0],param[1],param[2],\
                nucleus[nuc], x[nuc], data[nuc], x_GT[nuc], x_FF_1m[nuc], x_FF_0m[nuc], x_FF_2m[nuc],\
                 x_GT_EC[nuc], x_FF_1m_EC[nuc], x_FF_0m_EC[nuc], x_FF_2m_EC[nuc] ))
            else:
                file.write("{} {} {} {} {:e} {:e} {:e} {:e} {:e} {:e} {:e} {:e} {:e} {:e}\n".format(param[0],param[1],param[2],\
                nucleus[nuc], x[nuc], data[nuc], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))


    #sys.exit(-1)
    return


def get_me_response_K(param,nucleus,Jpi,k_val,types,recalc):
    '''
     calculates response of a particular nucleus

     Jpi is the angular momentum number defined above

     K is provided as an input in this module

     It is only used to correct for solutions that have not converged
     so I have commented part for the calculation of g.s.

    '''

    nucnam=nucleus[0:2]
    nucmas=nucleus[2:]

    # select isoscalar pairing strength
    V0_pp = param[0]
    g_A = param[1]
    g_0 = param[2]

    # calculate the ground-state
    # do it only for the first iteration

    # calculate the GT response
    calculate_response_K(nucleus,V0_pp,g_0,Jpi,k_val,types,recalc)




    return


# calculate GT resonance energies for spherical nuclei
def calculate_GT_resonance(param, nucleus):
    '''
      calculates GT resonance energies using the
      spherical_QRPA code

      param[0] - V_0^pp
      param[1] - g_A (not used)
      param[2] - g_0'

      nucleus - name of the nucleus

    '''

    # first run the ground-state calculations
    # to create the dir.wav file

    #Read the content of the file
    file_path='/mnt/gs21/scratch/ravlican/spherical_QRPA/code_FTPCRHB/dis.dat'
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Modify the specific line
    nucleus = nucleus.upper() # Pb208 --> PB208
    new_line_content = '{}                              ! nucleus under consideration\n'.format(nucleus)
    lines[7 - 1] = new_line_content  # Subtract 1 as Python uses 0-based indexing


    # Write the updated content back to the file
    with open(file_path, 'w') as file:
        file.writelines(lines)


    # change directory to GS code
    os.chdir('/mnt/gs21/scratch/ravlican/spherical_QRPA/code_FTPCRHB')
    subprocess.run(['./run_me.sh'], shell=True)

    # move dir.wav file to the QRPA folder
    # prepare parent directory
    folder_path='/mnt/gs21/scratch/ravlican/fit_me_beta_pm/chi2_fit/'
    os.makedirs(folder_path+nucleus, exist_ok=True)


    # now move the qrpa.wel file into the directory for that nucleus
    subprocess.run(['mv dir.wav ../../fit_me_beta_pm/chi2_fit/'+nucleus], shell=True)

    # now go to the subfolder with the nucleus
    os.chdir('/mnt/gs21/scratch/ravlican/fit_me_beta_pm/chi2_fit/'+nucleus)

    # copy files to run the sphrical QRPA code
    subprocess.run(['cp /mnt/gs21/scratch/ravlican/spherical_QRPA/code_FTQRPA_RHB/run .'], shell=True)
    subprocess.run(['cp /mnt/gs21/scratch/ravlican/spherical_QRPA/code_FTQRPA_RHB/rpa.dat .'], shell=True)
    subprocess.run(['cp /mnt/gs21/scratch/ravlican/spherical_QRPA/code_FTQRPA_RHB/rpa.par .'], shell=True)
    subprocess.run(['cp /mnt/gs21/scratch/ravlican/spherical_QRPA/code_FTQRPA_RHB/run_me.sh .'], shell=True)

    # modify the rpa.dat file with the new parameters
    file_path= 'rpa.dat'
    with open(file_path, 'r') as file:
        lines = file.readlines()

        # Modify the specific line
    new_line_content = 'pairfac  = {}                     ! isoscalar pairing strength\n'.format(round(param[0],5))
    lines[9 - 1] = new_line_content

    new_line_content = 'deltafac = {}                   ! landau-migdal parameter\n'.format(round(param[2],5))
    lines[10 - 1] = new_line_content

    # Write the updated content back to the file
    with open(file_path, 'w') as file:
        file.writelines(lines)

    # at this point you can sbatch the files
    subprocess.run([' ./run_me.sh'], shell=True)


    # in the updated version of the spherical code, I have moved the peak finidng
    # to F77 codes
    omega_cent = float(np.loadtxt('peak.dat'))

    #print('GTR energy found at:' , data[0][0], ' for ', nucleus, param, flush=True)
    print('GTR energy found at:' , omega_cent, ' for ', nucleus, param, flush=True)






    #return data[0][0]
    return omega_cent



def calculate_residuals_Jacobian(param, user_data):
    '''
    Helper function for calculate_Jacobian function

    input:
        param: arary of parameters
        user_data: data for calculation of residuals

    output:
        x: array of residuals

    '''

    global iteration_num

    data, nucleus, types = user_data # unpack tuple

    n_d = len(nucleus) # number of nuclei
    n_x = len(param) # number of parameters

    # I have to calculate all half-lives for a given
    # value of the parameter
    x = [] # list of calculated total half-lives
    E_GTR = [] # list of resonance energies
    for nuc in range(len(nucleus)):

        if (types[nuc] == 'HL-'):
            # gets to the queue
            for k in range(2):
                get_me_response_K(param,nucleus[nuc],'1+',k,'b-',False)
            for k in range(2):
                get_me_response_K(param,nucleus[nuc],'1-',k,'b-',False)
            
            get_me_response_K(param,nucleus[nuc],'0-',0,'b-',False)

            for k in range(3):
                get_me_response_K(param,nucleus[nuc],'2-',k,'b-',False)

        elif (types[nuc] == 'HL+'):

            for k in range(2):
                get_me_response_K(param,nucleus[nuc],'1+',k,'b+',False)
            for k in range(2):
                get_me_response_K(param,nucleus[nuc],'1-',k,'b+',False)
            
            get_me_response_K(param,nucleus[nuc],'0-',0,'b+',False)

            for k in range(3):
                get_me_response_K(param,nucleus[nuc],'2-',k,'b+',False)

            for k in range(2):
                get_me_response_K(param,nucleus[nuc],'1+',k,'EC',False)
            for k in range(2):
                get_me_response_K(param,nucleus[nuc],'1-',k,'EC',False)
            
            get_me_response_K(param,nucleus[nuc],'0-',0,'EC',False)

            for k in range(3):
                get_me_response_K(param,nucleus[nuc],'2-',k,'EC',False)


        else:
            # calculates straight away
            E_GTR.append(calculate_GT_resonance(param,nucleus[nuc]))

    #block untill all processes are finished
    block_jobs()

    # check that all jobs have converged before moving on
    check_jobs(param,user_data)





    print('All files have been checked and I am moving on to calculate HLs ...')


    # now calculate the half_lives
    for nuc in range(len(nucleus)):

        if (types[nuc] == 'HL-'):
            fac1 = GT_rate_module.calculate_HL_1p(nucleus[nuc],param[1])
            fac2 = FF_1m_rate_module.calculate_HL_1m(nucleus[nuc],param[1])
            fac3 = FF_0m_rate_module.calculate_HL_0m(nucleus[nuc],param[1])
            fac4 = FF_2m_rate_module.calculate_HL_2m(nucleus[nuc],param[1])
            fac5 = 0.0
            fac6 = 0.0
            fac7 = 0.0
            fac8 = 0.0
        elif (types[nuc] == 'HL+'):
            fac1 = GT_rate_module_Bp.calculate_HL_1p(nucleus[nuc],param[1])
            fac2 = FF_1m_rate_module_Bp.calculate_HL_1m(nucleus[nuc],param[1])
            fac3 = FF_0m_rate_module_Bp.calculate_HL_0m(nucleus[nuc],param[1])
            fac4 = FF_2m_rate_module_Bp.calculate_HL_2m(nucleus[nuc],param[1])
            fac5 = GT_rate_module_EC.calculate_HL_1p(nucleus[nuc],param[1])
            fac6 = FF_1m_rate_module_EC.calculate_HL_1m(nucleus[nuc],param[1])
            fac7 = FF_0m_rate_module_EC.calculate_HL_0m(nucleus[nuc],param[1])
            fac8 = FF_2m_rate_module_EC.calculate_HL_2m(nucleus[nuc],param[1])

        # now here correct for negative half-lives
        '''
            this can happen in the instances where V0_pp is too strong
            just assign 1e-10 to the rate in that case so that it is
            negligible

        '''

        if (fac1 <= 0.0):
            fac1 = 1e10
        if (fac2 <= 0.0):
            fac2 = 1e10
        if (fac3 <= 0.0):
            fac3 = 1e10
        if (fac4 <= 0.0):
            fac4 = 1e10
        if (fac5 <= 0.0):
            fac5 = 1e10
        if (fac6 <= 0.0):
            fac6 = 1e10
        if (fac7 <= 0.0):
            fac7 = 1e10
        if (fac8 <= 0.0):
            fac8 = 1e10


        # append the partial half-lives
        if (types[nuc] == 'HL-'):


            # append the total half-life
            x.append(1./(1./fac1 + 1./fac2 + 1./fac3 + 1./fac4))
        
        elif (types[nuc] == 'HL+'):


            # append the total half-life
            x.append(1./(1./fac1 + 1./fac2 + 1./fac3 + 1./fac4 + \
            1./fac5 + 1./fac6 + 1./fac7 + 1./fac8))

            


    # append the resonances at the end of the list
    for ee in range(len(E_GTR)):
        x.append(E_GTR[ee])




    # increase the number of iterations by 1
    iteration_num += 1
    print('Iteration number: ', iteration_num)



    return x


def calculate_Jacobian(optimal_param, user_data):
    '''
    Calculates the Jacobian matrix according to the formula:

    J_kl(x) = (res_k(x + d e_l) - res_k(x - d e_l))/2d

    where d = 1e-3, res_k is the residual, and e_l is the unit
    vector in the direction of the parameter.

    If the number of data points is n_d, and the number of parameters in n_x
    then the dimension of the Jacobian matrix is n_d x n_x

    input:
        optimal_param: array of optimal parameters from chi2 fit
        user_data: pass data for the residuals


    '''
    eps = 1e-3 # this is for finite-difference

    data, nucleus, types = user_data # unpack tuple

    # create array for Jacobian calculation
    for p in range(len(optimal_param)):

        param = np.array(optimal_param)
        param[p] = param[p] + eps
        param_plus = param

        param = np.array(optimal_param)
        param[p] = param[p] - eps
        param_minus = param

        print('Parameter check: ', p, optimal_param, param_plus, param_minus)
        #sys.exit(-1)


        x1 = np.array(calculate_residuals_Jacobian(param_plus, user_data))
        x2 = np.array(calculate_residuals_Jacobian(param_minus, user_data))

        # arrays of residuals
        res1 = np.zeros_like(x1)
        res2 = np.zeros_like(x2)


        # x1 and x1 contain HLs and GTR energies
        # I need to get residuals instead

        n_d = len(x1)
        n_x = len(optimal_param)

        for i in range(n_d):
        
            if (types[i] == 'HL-' or types[i] == 'HL+'):
                res1[i] = ((np.log10(x1[i])-np.log10(data[i]))/0.5)
                res2[i] = ((np.log10(x2[i])-np.log10(data[i]))/0.5)
            else:
                res1[i] = ((x1[i]-data[i])/0.30)
                res2[i] = ((x2[i]-data[i])/0.30)



        J = (res1 - res2)/2/eps


        # write the results of this Jacobian
        file_path = '/mnt/gs21/scratch/ravlican/fit_me_beta_pm/chi2_fit/Jacobian_'+str(p)+'.txt'

        # Open the file in write mode ('w')
        with open(file_path, 'w') as file:
            # write calculated HL for each nucleus
            file.write("# {} {} {} {}\n".format('plus', param_plus[0], param_plus[1], param_plus[2]))
            file.write("# {} {} {} {}\n".format('minus', param_minus[0], param_minus[1], param_minus[2]))

            for nuc in range(len(nucleus)):
                file.write("{} {:e} {:e} {:e}\n".format(nucleus[nuc], res1[nuc], res2[nuc], J[nuc]))



    return


if __name__ == "__main__":  

    iteration_num = 0



    # user selects nucleus for which to calculate HL
    nuclei = ['Ti58','Cr60','Zn78','Ce152',
    'Ar48','Kr96','Sn134','Gd166','Zr100',
    'Mo108','_U242','Xe142','Zn82','Ru114',
    'Nd156','Th236','Fe66','Pd120',
    'Sm162','Se88','Sr102','Rn228','Sr96','Er174','Pt204'
    ,'Fe50','Nd132','Yb160','Mo86','Zn60', 
    'Pb208','Sn132','Zr90','Sn112']

    t_exp=[0.058,0.49,1.47,1.4,0.416,
    0.080,1.050,4.8,7.1,1.09,1008,1.23,
    0.166,0.54,5.26,2250,0.440,0.492,
    2.4,1.53,0.069,65,1.07,192,10.3,
    0.152,94,288,19.1,143,
    19.2,13.2,15.5,16.8]

    # HL- means b- decay half-life
    # HL+ ,meanbs EC/b+ decay half-life

    types = ['HL-','HL-','HL-','HL-','HL-','HL-','HL-','HL-','HL-',
    'HL-','HL-','HL-','HL-','HL-','HL-','HL-','HL-','HL-','HL-',
    'HL-','HL-','HL-','HL-','HL-','HL-',
    'HL+','HL+','HL+','HL+','HL+'
    ,'GTR','GTR','GTR','GTR']


 


    # Create a TAO solver
    #tao = PETSc.TAO().create(PETSc.COMM_WORLD)

    ## set solver as POUNDERS
    #tao.setType(PETSc.TAO.Type.POUNDERS)

    ## Set the objective function
    #res = PETSc.Vec().createSeq(len(nuclei))
    #user_params = (t_exp,nuclei,types)
    #tao.setResidual(residual_routine,res,user_params)

    ##set variable bounds
    ## [V_0^pp, g_A, g_0']
    #lb = PETSc.Vec().createWithArray([0.0, 0.1,0.1])
    #ub = PETSc.Vec().createWithArray([1.8, 2.0,1.0])
    #tao.setVariableBounds([lb,ub])

    ## create initial guess vector
    #initial_guess = [1.0, 1.0, 0.734]
    ##initial_guess = [0.7540895130798597,0.6272078050735158,0.6293783091025107]
    #x = PETSc.Vec().createWithArray(initial_guess)
    #tao.setSolution(x)

    ## Optimize
    #tao.solve()

    #tao.view()

    ## Get the solution
    #solution = tao.getSolution()

    #print('---------------------------------', flush=True)
    #print("Optimal Solution found:", solution.getArray(), flush=True)
    #print('Number of iterations: ', iteration_num, flush=True)

    # calculate_GT_resonance([0.6864383790733146, 0.6681640625, 0.6565702755601917], 'Pb208')
    # sys.exit(-1)


    ## this is the part which just calculates the Jacobian
    user_params = (t_exp,nuclei,types)
    #opt_params = [1.41861281, 0.47423738, 0.47675899]
    opt_params = [1.431,0.704,0.613]
    calculate_Jacobian(opt_params, user_params)







     

    







