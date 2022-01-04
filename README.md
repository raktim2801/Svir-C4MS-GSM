# Svir-C4MS-GSM

_____________
#### Installation and requirements
We have used __scobra__ metabolic modelling package for solving our linear programming problems.
The specific version (usable in python 2.7) used is available as __supplementary material 8__ of the following [article](https://doi.org/10.3389/fpls.2018.00884) [[1]](#1).

###### Article
<a id='1'>[1]</a>
Shaw, R. and Cheung, C.Y., 2018. A dynamic multi-tissue flux balance model captures carbon and nitrogen metabolism and optimal resource partitioning during Arabidopsis growth. Frontiers in plant science, 9, p.884.

Latest version (works with python 3.6 onwards) of __scobra metabolic modelling__ package is available in [Github](https://github.com/mauriceccy/scobra).

This project utilizes __scobra__, __cobra__, __os__, __pickle__, __csv__, __numpy__ and __pandas__ modules.

_____________
_____________

#### Simulation
>All input files (Model (in excel format) and Gene Expression libraries are stored in
>__'.\\Data'__ folder)
>
>__Model File__:
>  - '.\\Data\\setaria_c4x4_v12a.xls'
>
>__Gene Expression values__:
>  -  '.\\Data\\GSMx2rxn_baseGE_N0_edited.txt',
>  -  '.\\Data\\GSMx2rxn_lwrGE_N0_edited.txt',
>  -  '.\\Data\\GSMx2rxn_midGE_N0_edited.txt',
>  -  '.\\Data\\GSMx2rxn_tipGE_N0_edited.txt'


__Simulation Guide__: (also given in **jupyter notebook** format) 'Model_simulation_guide.ipynb'


To start simply import the _setaria_final.py_ file.
```python
import setaria_final
```
it imports all necessary Modules on its own.

Load the model using..
```python
m=setaria_final.scobra.Model('.\\Data\\setaria_c4x4_v12a.xls')
```
Run the simulation using..

```python
MwGE, SolwGE = setaria_final.TwoStepBiomassGradientSolve(m, P=200.0, DWt='Mean',
                                                            r_tmlb=[1.283,1.2334,1.9372,5.488],
                                                            EqualPhoton=True)
```
By default it returns a model object and a scobra solution object.

The scobra solution object can be exported to csv file using..
```python
with open('wGE_Sol.csv','w') as outfile:
    writer = setaria_final.csv.writer(outfile, delimiter='\t')
    for key, value in SolwGE.items():
        writer.writerow([key, value])
```

If _SaveModel_, _SaveSols_ flags are __True__ in ```TwoStepBiomassGradientSolve```, it will create a __Results__ folder inside working directory and save (along with ```pickle``` dump) the models or the solutions or both inside that folder.

The ```pickle``` files (extension _.pkl_) can be loaded using the following commands..
```python
with open('.\\Results\\wGE_Model_run6.pkl', 'rb') as ff:
    wGESolved_M = setaria_final.pickle.load(ff)
```
