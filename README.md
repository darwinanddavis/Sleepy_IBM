<img src="https://raw.githubusercontent.com/darwinanddavis/Sleepy_IBM/master/img/sleepyibm_header.png" alt=" " width=1000 height=600>  

# Individual-based simulation model for space use and home range patterns under movement constraints based on location data of the model species *Tiliqua rugosa*.       

## Matthew Malishev<sup>1,2*</sup>, C. Michael Bull<sup>3</sup> & Michael R Kearney<sup>2</sup>   

#### _<sup>1</sup> Centre of Excellence for Biosecurity Risk, Melbourne, Victoria, Australia, 3010_      

#### _<sup>2</sup> School of Biological Sciences, University of Melbourne, Melbourne, Victoria, Australia, 3010_      

#### _<sup>3</sup> School of Biological Sciences, Flinders University, Adelaide, Australia, 5042_      

#### *Corresponding author: matthew.malishev [at] gmail.com     

## Overview  

Individual-based simulation model for forecasting home range patterns and dispersal potential of animals under food and microclimate constraints in space and time. The individual internal _i_ state is updated at each time step of the model with a dynamic energy budget (DEB) model to compute the current energy reserves, somatic maintenance costs, movement costs, growth in mass, and reproductive output potential. The microclimate of the habitat is computed at each time step from the [`NicheMapR`](https://github.com/mrke/NicheMapR) microclimate model based on hourly air temperature, air humidity, soil temperature, cloud coverage, and solar radiation data from the geolocation. The individual's energy, water, and mass balance updates according to changes in these combined internal and external cues to inform future space use.     
## Project outputs      

**Malishev M**, Michael Bull C, Kearney MR (2018) An individual‐based model of ectotherm movement integrating metabolic and microclimatic constraints. Methods Ecol Evol. 9:472–489. [https://doi.org/10.1111/2041-210X.12909](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12909#).    

[Supplementary material including model code and instructions.](https://github.com/darwinanddavis/MalishevBullKearney)       

*Files:*        

.R      
.nlogo  
.csv        

## Instructions  

* Download all files to your local drive.   
* Open the `RNL_new trans model_with DEB_1.6.2.R` in `R` and run the model code after setting your working directory.  

## Maintainer  
**Matt Malishev**   
:mag: [Website](https://www.researchgate.net/profile/Matt_Malishev)    
:bird: [@darwinanddavis](https://twitter.com/darwinanddavis)  
:email: matthew.malishev [at] gmail.com    
