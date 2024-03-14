
# To fill or not to fill: comparing imputation methods for improved riverine long-term biodiversity monitoring
## Leader: <a href="https://philliphaubrock.wixsite.com/invasivespecies">Phillip J. Haubrock</a>

<strong>Authors</strong>:
Phillip J. Haubrock1,2,3, Ismael Soto2, Rafael L. Macêdo4,5,6

<strong>Affiliations</strong>:

- **1** - Department of River Ecology and Conservation, Senckenberg Research Institute and Natural History Museum Frankfurt, Gelnhausen, Germany
- **2** - Faculty of Fisheries and Protection of Waters, South Bohemian Research Center of Aquaculture and Biodiversity of Hydrocenoses, University of South Bohemia in České Budějovice, Vodňany, Czech Republic
- **3** - CAMB, Center for Applied Mathematics and Bioinformatics, Gulf University for Science and Technology, Kuwait
- **4** - Institute of Biology, Freie Universität Berlin, Königin-Luise-Str. 1-3, 14195 Berlin, Germany
- **5** - Leibniz Institute of Freshwater Ecology and Inland Fisheries (IGB), Müggelseedamm 310, 12587 Berlin, Germany
- **6** - Graduate Program in Ecology and Natural Resources, and Department of Ecology and Evolutionary Biology, Federal University of São Carlos, Sao Carlos, Brazil


## Abstract
In the era of globalisation and climate change, the preservation of global aquatic biodiversity has become unprecedentedly challenging due to intensifying and increasingly synergistic anthropogenic pressures. This study addresses the complex challenges associated with long-term biodiversity monitoring data, specifically the effects of missing years and often applied gap-filling on the accuracy of inferring temporal trends of riverine macroinvertebrate biodiversity. For this, we utilised three distinct time series of >20 years of annual riverine macroinvertebrate community samplings from Denmark, the Netherlands, and Sweden. We apply random deletions to simulate incomplete data and employ multiple imputation methods, including Predictive Mean Matching, Weighted Predictive Mean Matching, Random Forest Imputations, and Random Sample from Observed Values, to fill gaps using linear and non-linear models. Our findings indicate that as the number of gaps increases, variability in trends increases concomitant to a decrease in the Akaike Information Criterion and rise in the standard deviation of the explained deviance by models. This indicates easier fitting with increasing data gaps and an elevated level of unexplained variability introduced by these gaps. Assessing the performance of gap-filling algorithms, we observed that Predictive Mean Matching, which is the most commonly employed imputation algorithm, exhibits limitations, leading to increased uncertainty. Random Forest Imputations and Random Sample from Observed Values demonstrated resilience in capturing underlying linear and non-linear temporal patterns, outperforming Predictive Mean Matching. The study underscores the nuanced challenges associated with missing data and its interpretation, emphasising the importance of careful modelling approach in biodiversity trend analyses. 

**Keywords**: Biodiversity monitoring; Missing data; Gap filling; Time series analysis; Resilience; Global biodiversity trends

## <a href="https://github.com/IsmaSA/Aquaculture/tree/master/Code">Script:</a>
Code created by: <a href="https://philliphaubrock.wixsite.com/invasivespecies">Phillip J. Haubrock</a>

Reproducible to create the plots of the paper
- <code>Fran.R</code>
