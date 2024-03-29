# SVJL helper

Package to generate self-consistent configuration cards for the pythia 8 hidden valley module for [SVJL signatures](https://link.springer.com/article/10.1140/epjc/s10052-022-10775-2). Please cite this paper when using this package, in addition to the theory papers included in the "citation" string output of the code and the repository citation (in the citation button and citation output string).

The code is an adapted version of: [CMS SVJ production](https://github.com/cms-svj/SVJProduction) , [dark_showers_tool](https://gitlab.com/simonknapen/dark_showers_tool).

## Usage
The code takes the following arguments:
  * ```--mZprime```: Z' mass (GeV) (required)
  * ```--rinv```: invisible fraction (required)
  * ```--svjl_type```: Versions of SVJL model (default: A-Democratic, as in  [SVJL signatures](https://link.springer.com/article/10.1140/epjc/s10052-022-10775-2)) (optional)
  * ```--mPiOverLambda```: lightest dark pseudoscalar mass over LambdaHV (required)
  * ```--lambda```: dark sector confinement scale (GeV) (required)
  * ```--nevents```: Number of events to generate (required)
  * ```--card_author```: author of the generated card (optional)
    
An example card is provided in the repository, in order to produce your own card use commands like:

```python svjHelper.py --mZprime 3000.0 --rinv 0.3 --mPiOverLambda 1.6 --lambda 5 --nevents 1000```


## Citations

The code itself can be cited as

```
 @software{Cazzaniga_SVJL_models_production_2023,
            author = {Cazzaniga, Cesare},
            month = dec,
            title = {{SVJL models production}},
            url = {https://github.com/cesarecazzaniga/svjl_helper},
            version = {1},
            year = {2023}
         }
```

## Support

If you have any questions, please
write us at 

[cesare.cazzaniga@cern.ch](cesare.cazzaniga@cern.ch)
