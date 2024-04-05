# 3D-CMCC-FEM  
# *Three Dimensional - Coupled Model Carbon Cycle - Forest Ecosystem Module*

![FML_4](https://user-images.githubusercontent.com/27204956/99875955-c8bfea80-2bf3-11eb-821c-c02a068990d4.jpg)

**3D-CMCC-FEM *(Three Dimensional - Coupled Model Carbon Cycle - Forest Ecosystem Module)* official repository**

LICENSE:

Copyright(c) 2020, Forest Modelling Laboratory – 3D-CMCC-FEM  
All rights reserved


The **3D-CMCC-FEM** is basically a research tool which is freely available ***only*** for non-commercial use. We have developed the **3D-CMCC-FEM** code relying solely on open source components, in order to facilitate its use and further development by others. 

The **3D-CMCC-FEM** is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. The **3D-CMCC-FEM** code is released under the GNU General Public Licence v3.0 (GPL). See the LICENSE file along with this program for more details.  If not, see http://www.gnu.org/licenses/gpl.html.

This page contains all the code releases developed over the time at CNR on the open source distribution of the computer simulation forest model **3D-CMCC-FEM**. The model has been developed and is maintained by the **Forest Modelling Laboratory** at the National Research Council of Italy, Institute for Agricultural and Forestry Systems in the Mediterranean (CNR-ISAFOM), Perugia (Italy). 

All source code and documents provided here are subject to copyright(c) by the CNR

In case you have copied and/or modified the **3D-CMCC-FEM** code overall, even in small parts of it, you may not publish data from it using the name **3D-CMCC-FEM** or any **3D-CMCC-FEM** variants unless you have either coordinated your usage and their changes with the developers listed below, or publish enough details about your changes so that they could be replicated.

The **3D-CMCC-FEM** has been developed by:  
Alessio Collalti  
Daniela Dalmonech  
Alessio Ribeca     
 


*Forest Ecology Laboratory* at the National Research Council of Italy (CNR),Institute for Agricultural and Forestry Systems in the Mediterranean (ISAFOM), Via della Madonna Alta, 128, 06128 - Perugia (PG), Italy.   

*DISCLAIMER*  
CNR accepts no responsibility for the use of the **3D-CMCC-FEM** in the form supplied or as subsequently modified by third parties. CNR disclaims liability for all losses,damages and costs incurred by any person as a result of relying on this software. Use of this software assumes agreement to this condition of use. Removal of this statement violates the spirit in which **3D-CMCC-FEM** was released by CNR.

For any further information please contact us:  
-Alessio Collalti (Lab. Head and Model Principal Investigator): alessio.collalti@cnr.it  
-Daniela Dalmonech: daniela.dalmonech@cnr.it  
-Mauro Morichetti: mauro.morichetti@cnr.it  
-Elia Vangi: elia.vangi@isafom.cnr.it  
-Elisa Grieco: elisa.grieco@cnr.it  
-The model: 3d-cmcc-fem@gmail.com  
-The lab: forest.modelling.lab@gmail.com  

or look at the **3D-CMCC-FEM** web-site:https://www.forest-modelling-lab.com/

For references on **3D-CMCC-FEM**:

1. *"Sviluppo di un modello ecologico-forestale per foreste a struttura complessa"*. Collalti A., Valentini R.. https://dspace.unitus.it/handle/2067/2398, *PhD Thesis*, University of Tuscia, 2011

2. *“A process-based model to simulate growth in forests with complex structure: Evaluation and use of 3D-CMCC Forest Ecosystem Model in a deciduous forest in Central Italy“*. Collalti A., Perugini L., Santini M., Chiti T., Nolè A., Matteucci G., Valentini R., *Ecological Modelling*, 272, 362– 378, https://doi.org/10.1016/j.ecolmodel.2013.09.016, 2014.

3. *“Assessing NEE and Carbon Dynamics among 5 European Forest types: Development and Validation of a new Phenology and Soil Carbon routines within the process oriented 3D-CMCC-Forest-Ecosystem Model“*. Marconi S., Chiti T., Valentini R., Collalti A., *MsC Master Thesis (U.S. – Italy Fulbright programme winner)*, University of Tuscia,  2014.

4. *”Validation of 3D-CMCC Forest Ecosystem Model (v.5.1) against eddy covariance data for ten European forest sites”*, A. Collalti, S. Marconi, A. Ibrom, C. Trotta, A. Anav, E. D'Andrea, G. Matteucci, L. Montagnani, B. Gielen, I. Mammarella, T. Grünwald, A. Knohl, F. Berninger, Y. Zhao, R. Valentini, and M. Santini. *Geoscientific Model Development*, 9:479–504, https://doi.org/10.5194/gmd-9-479-2016, 2016.

5.  *“Protocollo di simulazione, calibrazione e validazione del modello 3D-CMCC-CNR FEM: il caso studio del bacino altamente strumentato del Bonis in Calabria”*, A. Collalti,  Biondo C., Buttafuoco G., Maesano M., Caloiero T., Lucà F., Pellicone G., Ricca N., Salvati R., Veltri A., Scarascia Mugnozza. G., Matteucci G.; *Forest@*, 14:247-256, https://doi.org/10.3832/efor2368-014, 2017.

6. *“The role of Respiration in estimation of the net Carbon cycle: coupling soil Carbon dynamics and canopy turnover in a novel version of 3D-CMCC Forest Ecosystem Model”*, Marconi S., Chiti T., Valentini R., Collalti A.. *Forests*, 8:220, https://doi.org/10.3390/f8060220, 2017.

7. *“Climate change mitigation by forests: a case study on the role of management on carbon dynamics of a pine forest in South Italy”* Pellicone G., Scarascia-Mugnozza G., Matteucci G., Collalti. A. https://dspace.unitus.it/handle/2067/42735, *PhD Thesis*, University of Tuscia, 2018.

8. *“Thinning can reduce losses in carbon use efficiency and carbon stocks in managed forests under warmer climate”*. A. Collalti , C. Trotta, T. Keenan, A. Ibrom, B. Bond-Lamberty, B. R. Grote, S. Vicca, C.P.O. Reyer, M. Migliavacca, F. Veroustraete, A. Anav, M. Campioli, E. Scoccimarro, E. Grieco, A. Cescatti, G. Matteucci. *Journal of Advances in  Modelling Earth Systems*, 10(10):2427-2452, https://doi.org/10.1029/2018MS001275, 2018.

9. *“The sensitivity of the forest carbon budget shifts between different parameters and processes along stand development and climate change”*, Collalti A., Thornton P. E., Cescatti A., Rita A., Nolè A., Borghetti M., Trotta C., Ciais P., Matteucci G..*Ecological Applications*, 29(2):1-18, https://doi.org/10.1002/eap.1837, 2019.

10. *“Forest carbon allocation modelling under climate change”*, Merganičová K., Merganič J., Lehtonen A., Vacchiano G., Ostrogović Sever M.Z., Augustynczik A.L.D., Grote R., Kyselova I., Mäkelä A., Yousefpour R., Krejza J., Collalti A., Reyer C.P.O.. *Tree Physiology*, 39:1937–1960, https://doi.org/10.1093/treephys/tpz105, 2019.

11. *“Plant respiration: Controlled by photosynthesis or biomass?”*, Collalti A., Tjoelker M.G., Hoch G., Mäkelä A., Guidolotti G., Heskel M., Petit G., Ryan M.G., Battipaglia G., Matteucci G., Prentice I.C.. *Global Change Biology*, 26: 1739–1753, https://doi.org/10.1111/gcb.14857, 2020.

12. *“Simulating the effects of thinning and species mixing on stands of oak (Quercus petrea (Matt.) Liebl. / Quercus robur L.) and pine (Pinus sylvestris L.) across Europe"*, Vorspernik E.M., Toigo M., Morin X., Tomao A., Trotta C., Steckel M., Barbati A., Nothdurft A., Trezsch H, del Rio M., Skrzyszewski J., Ponette Q., Lof M., Jansons A., Brazaitis G. *Ecological Modelling*, https://doi.org/10.1016/j.ecolmodel.2020.109406, 2021.
    
13. *“Accuracy, realism and general applicability of European forest models"*, Mahnken M., Cailleret M., Collalti A., Trotta C., Biondo C., D'Andrea E., Dalmonech D., Marano G., Mäkelä A., Minunno F., Peltoniemi M., Trotsiuk V., Nadal-Sala D., Sabaté S., Vallet P., Aussenac R., Cameron D.R., Bohn F.J., Grote R., Augustynczik A.L.D., Yousefpour R., Huber N., Bugmann H., Merganikova K., Merganic J., Valent P., Lash-Born P., Hartig F., Vega del Valle I.D., Volkholz J., Gustch M., Matteucci G., Krejza J., Ibrom A., Meesenburg H., Rötzer T., van der Maaten-Theunissen M., van der Maaten E., Reyer C.P.O. *Global Change Biology*,  https://doi.org/10.1111/gcb.16384, 2022.

14. *“Feasibility of enhancing carbon sequestration and stock capacity in temperate and boreal European forests via changes to forest management"*, Dalmonech D. , Marano G. , Amthor J., Cescatti A., Lindner M., Trotta C., Collalti A.. *Agricultural and Forest Meteorology*, 327: 109203, https://doi.org/10.1016/j.agrformet.2022.109203, 2022.

15. *"3D-CMCC-FEM User’s Guide"*, Collalti A., Dalmonech D., Marano G., Vangi E., Puchi P., Grieco E., Orrico M., *CNR Edizioni*,  ISBN 978-88-8080-573-1, https://doi.org/10.32018/3D-CMCC-FEM-2022, 2022.

16. *“Simulating diverse forest management in a changing climate on a Pinus nigra subsp. Laricio plantation in Southern Italy"*, Testolin R., Dalmonech D., Marano G., D’Andrea E., Matteucci G., Noce S., Collalti A.. *Science of the Total Environment*, 857: 159361, https://doi.org/10.1016/j.scitotenv.2022.159361, 2023.

17. *“Regional estimates of gross primary production applying the Process-Based Model 3D-CMCC-FEM vs. Remote-Sensing multiple datasets"*, Dalmonech D., Vangi E., Chiesi M., Chirici G., Fibbi L., Giannetti F., Marano G., Massari C., Nolè A., Xiao J., Collalti A.. *European Journal of Remote Sensing*, 57(1): 2301657, https://doi.org/10.1080/22797254.2023.2301657, 2024.

18. *"A harmonized database of European forest simulations under climate change"*, Grünig M., Rammer W., Albrich A., Andrè F., Augustynczik A.L.D., Bohn F., Bowman M., Bugmann H., Collalti A., Cristal I., Dalmonech D., De Caceres M., De Coligny F., Dobor L., Dollinger C., Forrester D.I., Garcia-Gonzalo J., González J.R., Hiltner U., Hlásny T., Honkaniemi J., Huber N., Jonard M., Jönsson A.M., Lagergren F., Nieberg M., Mina M., Mohren F., Moos C., Morin X., Muys B., Peltoniemi M., Reyer C.P.O., Storms I., Thom D., Toigo M., Seidl R.. *Data in Brief*, 110384, https://doi.org/10.1016/j.dib.2024.110384, 2024.

see also at: https://www.forest-modelling-lab.com/publications

