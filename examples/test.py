import numpy

import pycrobe.objects as pycrobe




strains = {}
strains['WT'] = pycrobe.Strain( name=strainName, 
						max_growth_rate=growthRates[strainName], 
						mean_lag_exit_time=1.5,
						stdev_lag_exit_time=0.1,
						marker='gfp')
strains['M1'] = pycrobe.Strain( name=strainName, 
						max_growth_rate=growthRates[strainName], 
						mean_lag_exit_time=1.5,
						stdev_lag_exit_time=0.1,
						marker='gfp' )

strains