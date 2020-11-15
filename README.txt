Paper: On Null Space-Based Inverse Kinematics Techniques for Fleet Management: Toward Time-Varying
Task Activation

Authors: Anna Mannucci, Danilo Caporale, Lucia Pallottino.

Journal: T-RO

Year: 2020


----------------------------------
         FILE DESCRIPTION
----------------------------------

main.m                  The main file. It launches all the simulations with
                        different controllers and plots results. Run this file section by section
                        (all the errors are plot even if the related task is not activated).

entrapment.m            This file is called in the main. It initializes references, 
                        simulation parameters and scenario
                 
entrapment_multi.slx    Simulink file to run simulations. The core of the 
                        file is the NSB_IK_controller which should be customized 
                        by the used according to the set of tasks to be enabled.

UTILITIES

damped_pinv.m           Function to compute the damped pseudo inverse (see Section III.A)

pinvSC.m                Function to compute the regularized pseudo inverse (see Section IV).

activation_simple.m     Function to compute a smooth activation value.

backtr_circ.m,          Eq. 14-15 of the paper.
costfun_circular.m           

compute_centroid.m      Function to compute the centroid of the fleet.

handle_name.m           Mapping between controller names (string) and RPTflag in the slx file.


----------------------------------
   HOW TO CUSTOMIZE SIMULATIONS
----------------------------------
CONTROLLERS Change the variable tests in the main.m (line 22).

TASKS Changing the NSB_IK_controller in the entrapment_multi.slx as follows to define
the tasks you want to enable:

1) line 4 (number of tasks) -- n_tasks
2) line 13 (number of subtasks) -- h_length_
3) lines 92-104 to use the gradient for defining the circular task (see
   Section IV.C)
4) lines 143-176 (piling up the correct Jacobians)

PARAMETERS Damping parameters, references and scenario: customize them inside the entrapment.m file