import numpy
import pandas

import pycrobe.objects as pycrobe




strains = {}
strains['WT'] = pycrobe.Strain( name='WT', 
						max_growth_rate=1.32, 		# randomize
						mean_lag_exit_time=1.5,		# randomize
						stdev_lag_exit_time=0.1,	# randomize
						marker='gfp')
strains['MX'] = pycrobe.Strain( name='MX', 
						max_growth_rate=1.30, 		# randomize
						mean_lag_exit_time=1.5,		# randomize
						stdev_lag_exit_time=0.1,	# randomize
						marker='rfp' )

flaskMediaVolume = 25.0
competitionFlask = pycrobe.CultureVolume( media=pycrobe.Media(volume=flaskMediaVolume, carrying_capacity=numpy.random.normal(3.0e9, 1.0e9)),
										  inoculums=[ pycrobe.Inoculum(strain=strains['WT'], cell_count=1.0e5*flaskMediaVolume), 
													  pycrobe.Inoculum(strain=strains['MX'], cell_count=1.0e5*flaskMediaVolume) ] )

incubator = pycrobe.Incubator(set_temp=37.0, temp_std_batch=0.5, temp_std_location=0.25, temp_std_transient=0.01)

allCompetitionData_alleleFreqs = pandas.DataFrame(columns=['competition', 'competition_rep', 'time', 'sample_rep', 'allele', 'allele_freq'])
allCompetitionData_popSizes    = pandas.DataFrame(columns=['competition', 'competition_rep', 'time', 'sample_rep', 'population', 'pop_size', 'density'])
allCompetitionData_growthCurves = pandas.DataFrame(columns=['culture', 'competition_rep', 'time', 'temp', 'density']) 

# incubator.incubate(cultures=[competitionFlask], time=24, dt=0.01, record_growth_curves=True)

competitionFlask.info()

flow = pycrobe.FlowCytometer()



def collectEnumerationSampleData(competition_culture, sample_reps, t, transfer_dfs, data, sim_flow_cytometry=True):
    for r in range(sample_reps):
        enumerationSample         = competition_culture.sample(1.0)
        # enumerationSample.in fo()
        enumerationDilutionSeries = DilutionSeries(enumerationSample)
        
        flowReadVolume    = 0.100
        flowReadDilutions = [-2, -1]
        for df in flowReadDilutions:
            
            flowEvents = flow.read(sample=enumerationDilutionSeries.getDilution(df), read_volume=flowReadVolume)

            totalDensityEstimate = 0.0
            totalPopSizeEstimate = 0.0
            for markerName, numEvents in flowEvents.iteritems():
                markerDensityEstimate = numEvents*(1.0/flowReadVolume) * 10**(numpy.abs(df)) 
                markerPopSizeEstimate = markerDensityEstimate * 10**(numpy.sum(transfer_dfs))  # reduce((lambda x, y: x * y), [competitionFlaskMediaVolume/xfrvol for xfrvol in transfer_volumes], 1)  # (competitionFlaskMediaVolume/competitionFlaskDefaultTransferVolume)**transferCount
                totalDensityEstimate += markerDensityEstimate
                totalPopSizeEstimate += markerPopSizeEstimate
                data['pop_sizes'].append( {'time': t, 'sample_rep': r+1, 'dilution': df,
                                           'population': markerName, 'pop_size': markerPopSizeEstimate, 'density':markerDensityEstimate} )
            data['pop_sizes'].append( {'time': t, 'sample_rep': r+1, 'dilution': df,
                                       'population': 'total', 'pop_size': totalPopSizeEstimate, 'density':totalDensityEstimate} )
           


def collectSequencingSampleData(competition_culture, sample_reps, t, data, sequencing_error=0.01, frequency_precision=0.001):
    for r in range(sample_reps):
        sequencingSample         = competition_culture.sample( competition_culture.totalVolume() )
        
        # In lieu of simulating sequencing more directly, return estimated frequencies based on cell counts with small noise:
        trueCellCounts = sequencingSample.getCellCounts()
        estimatedCellCounts = {}
        for inoculumName, cellCount in trueCellCounts.iteritems():
            estimatedCellCounts[inoculumName] = int(numpy.random.normal(cellCount, cellCount*sequencing_error))
        estimatedTotalCellCount = numpy.sum(estimatedCellCounts.values())

        for inoculumName, cellCount in estimatedCellCounts.iteritems():
            frequencyEstimate = round( cellCount/estimatedTotalCellCount, int(numpy.log10(1/frequency_precision)) )
            data['allele_freqs'].append( {'time': t, 'sample_rep': r+1, 
                                          'allele': inoculumName, 'allele_freq': frequencyEstimate} ) 

competitionFlaskMediaVolume    = 25.0

totalTime            = 36
transferTimeInterval = 4

numCompetitionReps = 1
competitionRep = 1 #placeholder because this notebook is not set up to run multiple reps

t = 0

while t <= totalTime:
    print "\t\tT=" + str(t) + "h"     

    competitionTimepointData = {'allele_freqs':[], 'pop_sizes':[]}

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Sample Competition Culture for Population Size Estimation (via flow cytometry enumeration):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print "\t\t\t\tPop size sampling"
    collectEnumerationSampleData(competition_culture=competitionFlask, 
                                 sample_reps=3, 
                                 t=t, 
                                 transfer_dfs=transferDfs[competition][:-1],
                                 data=competitionTimepointData,
                                 sim_flow_cytometry=False)
    # Annotate this time point data with the competition name and rep number:
    competitionTimepointData['pop_sizes']    = [dict(d, competition="competition", competition_rep=1) for d in competitionTimepointData['pop_sizes']]
    # Add this competition's data recorded for this time point to the overall dataset:
    allCompetitionData_popSizes    = allCompetitionData_popSizes.append(competitionTimepointData['pop_sizes'])


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Sample Competition for Allele Frequency Estimation (via sequencing):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print "\t\t\t\tAllele freq sampling"
    collectSequencingSampleData(competition_culture=competitionFlask, 
                                sample_reps=1, 
                                t=t, 
                                data=competitionTimepointData)
    # Annotate this time point data with the competition name and rep number:
    competitionTimepointData['allele_freqs'] = [dict(d, competition="competition", competition_rep=1) for d in competitionTimepointData['allele_freqs']]
    # Add this competition's data recorded for this time point to the overall dataset:
    allCompetitionData_alleleFreqs = allCompetitionData_alleleFreqs.append(competitionTimepointData['allele_freqs'])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Incubate the Competition Flasks:
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print "\t\t\tIncubating flasks..."
    growthCurves = incubator.incubate(cultures=competitionFlask, time=transferTimeInterval, dt=0.01, record_growth_curves=True)
    growthCurves['time'] += t
    growthCurves['competition_rep'] = competitionRep
    allCompetitionData_growthCurves = allCompetitionData_growthCurves.append( growthCurves )
    

    if(t != 0):
        transferCount += 1


    # Update the time point
    t = t+transferTimeInterval