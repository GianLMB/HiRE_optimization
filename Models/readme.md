

Here ML models are reported, with corresponding datasets, source code, results and analyses. In particular, [Model_first](Model_first) contains the first implementation, with results obtained with classical SGD, Preconditioned Stochastic Gradient descent (PSGD) and fixation of geometric parameters; [Model_less_dof](Model_less_dof) eliminates degeneracy in some global parameters and [Model_with_tanh](Model_with_tanh) introduces different functional forms to impose soft constraints in geometric parameters.

All directories are structures in a similar way:
```
Model_*
│   
└───analyses 
│   |   jupyter notebooks
|
└───data
│   │   datasets
│   
└───results  
│   |   results.pth files  
│       
└───script      
│   |   script.py to execute
│   
└───src
    │   data_classes.py
    │   energy_functions.py
    |   optim_functions.py
    |   ...  all Python source code
```
In order to perform an optimization of chosen parameters, it is necessary to access the file script/script.py to modify the code. Results are then saved in results/ and can be analysed with the scripts contained in analyses.
