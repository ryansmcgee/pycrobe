from __future__ import division
from decimal import Decimal

import random
import numpy
import pandas
import scipy.stats
import copy


##################################################
##################################################


class Solute(object):

	def __init__(self, name="Solute", concentration=0, molecules_per_unit_conc=1e8):

		self.name 			= name

		self.molecules_per_unit_conc = molecules_per_unit_conc

		#------------------------------
		# Variables:
		#------------------------------
		self.concentration 	= concentration


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def sample(self, sample_volume, original_volume):
		sample_proportion 		= sample_volume / original_volume
		molecules_per_unit_vol 	= self.concentration * self.molecules_per_unit_conc 	# ug/mL * molecules/ug = molecules/mL
		originalNumMolecules 	= original_volume * molecules_per_unit_vol 		# mL * molecules/mL = molecules
		sampleNumMolecules 		= numpy.random.poisson(originalNumMolecules * sample_proportion)	# molecules

		sampleConcentration 	= (sampleNumMolecules/self.molecules_per_unit_conc) / sample_volume 	# (molecules / molecules/ug) / mL = ug/mL
		
		self.concentration 		= ((originalNumMolecules-sampleNumMolecules)/self.molecules_per_unit_conc) / (original_volume - sample_volume)	# (molecules / molecules/ug) / mL = ug/mL

		sampledSolute 				= copy.deepcopy(self)
		sampledSolute.concentration = sampleConcentration 

		return sampledSolute


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def add(self, added_solute, pct_of_current_volume_added):
		assert(isinstance(added_solute, type(self)) and self == added_solute), "To add to a Solute, the added object must have matching type and parameters."
		self.concentration 	= (self.concentration + added_solute.concentration*pct_of_current_volume_added) / (1 + pct_of_current_volume_added)
		return


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def __eq__(self, other):
		return(	self.name == other.name and type(Solute) == type(Solute))


##################################################
##################################################


class Nutrient(Solute):

	def __init__(self, name="Nutrient", concentration=0, molecules_per_unit_conc=1e8):

		super(Nutrient, self).__init__(name, concentration, molecules_per_unit_conc)


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	# Does superclass (Solute) sample method suffice???

	# def sample(self, sample_volume, original_volume):
	# 	sample_proportion 		= sample_volume / original_volume
	# 	molecules_per_unit_vol 	= self.concentration * molecules_per_unit_conc 	# ug/mL * molecules/ug = molecules/mL
	# 	originalNumMolecules 	= original_volume * molecules_per_unit_vol 		# mL * molecules/mL = molecules
	# 	sampleNumMolecules 		= numpy.random.poisson(originalNumMolecules * sample_proportion)	# molecules
		
	# 	sampleConcentration 	= (sampleNumMolecules/molecules_per_unit_conc) / sample_volume 	# (molecules / molecules/ug) / mL = ug/mL
	# 	self.concentration 		= ((originalNumMolecules-sampleNumMolecules)/molecules_per_unit_conc) / sample_volume 	# (molecules / molecules/ug) / mL = ug/mL

	# 	sampledNutrient 				= copy.deepcopy(self)
	# 	sampledNutrient.concentration 	= sampleConcentration 

	# 	return sampledNutrient


##################################################
##################################################


class Drug(Solute): 

	def __init__(self, name="Drug", concentration=0, decay_rate=0, molecules_per_unit_conc=1e8):

		super(Drug, self).__init__(name, concentration, molecules_per_unit_conc)

		#------------------------------
		# Parameters:
		#------------------------------
		self.decayRate		= decay_rate									# delta_A	


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	# Does superclass (Solute) sample method suffice???

	# def sample(self, sample_volume, original_volume):
	# 	sample_proportion 		= sample_volume / original_volume
	# 	molecules_per_unit_vol 	= self.concentration * molecules_per_unit_conc 	# ug/mL * molecules/ug = molecules/mL
	# 	originalNumMolecules 	= original_volume * molecules_per_unit_vol 		# mL * molecules/mL = molecules
	# 	sampleNumMolecules 		= numpy.random.poisson(originalNumMolecules * sample_proportion)	# molecules
		
	# 	sampleConcentration 	= (sampleNumMolecules/molecules_per_unit_conc) / sample_volume 	# (molecules / molecules/ug) / mL = ug/mL
	# 	self.concentration 		= ((originalNumMolecules-sampleNumMolecules)/molecules_per_unit_conc) / sample_volume 	# (molecules / molecules/ug) / mL = ug/mL

	# 	sampledDrug 				= copy.deepcopy(self)
	# 	sampledDrug.concentration 	= sampleConcentration 

	# 	return sampledDrug


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def __eq__(self, other):
		return(	self.name == other.name and type(self) == type(other) and self.decayRate == other.decayRate)


##################################################
##################################################


class Media(object):

	def __init__(self, volume, drugs=[], solutes=[], nutrient=Nutrient(), is_intracellular=False):

		self.volume		= volume

		self.drugs		= []	# class instance list has to be declared this way to avoid list being shared among instances when mutable changes are made to one instance list, apparently
		for drug in (drugs if isinstance(drugs, list) else [drugs]):
			self.drugs.append(drug)

		self.solutes	= []	# class instance list has to be declared this way to avoid list being shared among instances when mutable changes are made to one instance list, apparently
		for solute in (solutes if isinstance(solutes, list) else [solutes]):
			self.solutes.append(solute)					

		self.nutrient	= nutrient

		self.isIntracellular = is_intracellular


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def sample(self, sample_volume):
		
		sampleProportion 	= sample_volume / self.volume

		sampledDrugs 		= [ drug.sample(sample_volume=sample_volume, original_volume=self.volume) for drug in self.drugs ]
		sampledSolutes 		= [ solute.sample(sample_volume, self.volume) for solute in self.solutes ]
		sampledNutrient 	= self.nutrient.sample(sample_volume, self.volume)

		sampledMedia 		= Media(volume=sample_volume, drugs=sampledDrugs, solutes=sampledSolutes, nutrient=sampledNutrient, is_intracellular=self.isIntracellular)

		self.volume 		-= sample_volume

		return sampledMedia


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def add(self, added_media):
		currentVolume 		= self.volume
		addedVolume 		= added_media.volume
		addedVolumeRatio 	= addedVolume / currentVolume
		# Add the other media's volume to this volume:
		self.volume += addedVolume
		# Add the other media's solutes to this media:
		for addedDrug in added_media.drugs:
			matchingExistingDrug = next((drug for drug in self.drugs if drug == addedDrug), None)
			if(matchingExistingDrug):
				# A matching drug already exists in this media, add the two drugs together:
				matchingExistingDrug.add(addedDrug, addedVolumeRatio)
			else:
				# The added drug does not exist in this media, introduce a new drug object to the culture:
				newDrug = copy.deepcopy(addedDrug)
				newDrug.concentration=0.0
				newDrug.add(addedDrug, addedVolumeRatio)
				self.drugs.append( newDrug )
		#...
		for addedSolute in added_media.solutes:
			matchingExistingSolute = next((solute for solute in self.solutes if solute == addedSolute), None)
			if(matchingExistingSolute):
				# A matching solute already exists in this media, add the two solutes together:
				matchingExistingSolute.add(addedSolute, addedVolumeRatio)
			else:
				# The added solute does not exist in this media, introduce a new solute object to the culture:
				newSolute = copy.deepcopy(addedSolute)
				newSolute.concentration=0.0
				newSolute.add(addedSolute, addedVolumeRatio)
				self.solutes.append( newSolute )
		#...
		self.nutrient.add(added_media.nutrient, addedVolumeRatio)
		# Note: Added medias do not change the isIntracellular of this media.
		return


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def info(self):
		if(self.nutrient is not None):
			print "\tNutrient " + self.nutrient.name + ":\n" + "\t\tconcentration\t= " + ("%.2E" % self.nutrient.concentration) + " ug/mL"
		else:
			print "\t(no nutrient)"
		#----------
		if(len(self.drugs)>0):
			for drug in self.drugs:
				print "\tDrug " + drug.name + ":\n" + "\t\tconcentration\t= " + ("%.3f" % drug.concentration) + " ug/mL"
		else:
				print "\t(no drugs)"
		#----------
		if(len(self.solutes)>0):
			for solute in self.solutes:
				print "\tSolute " + solute.name + ":\n" + "\t\tconcentration\t= " + (("%.3f" % solute.concentration) if solute.concentration>0.001 else ("%.2E" % solute.concentration)) + " ug/mL"
		else:
				print "\t(no other solutes)"



##################################################
##################################################


class Strain(object):

	def __init__(self,name, max_growth_rate, optimal_temp=37.0, mean_lag_exit_time=1.5, stdev_lag_exit_time=0.125, 
							halfmax_growth_nutrient_conc=numpy.inf, nutrient_consumption_rate=0,
							marker='', plasmids=[]):

		self.name 								= name

		#------------------------------
		# Parameters:
		#------------------------------
		self.maxGrowthRate 						= max_growth_rate				# rho
		self.optimalTemp						= optimal_temp
		self.meanLagExitTime 					= mean_lag_exit_time
		self.stdevLagExitTime 					= stdev_lag_exit_time

		self.halfmaxGrowthNutrientConc 			= halfmax_growth_nutrient_conc 	# k_r		
		self.nutrientConsumptionRate 			= nutrient_consumption_rate		# sigma

		self.marker 							= marker

		self.plasmids 							= plasmids


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def __eq__(self, other):
		return(self.name == other.name 
				and type(self) 						== type(other)
				and self.maxGrowthRate 				== other.maxGrowthRate
				and self.optimalTemp				== other.optimalTemp
				and self.meanLagExitTime 			== other.meanLagExitTime
				and self.stdevLagExitTime 			== other.stdevLagExitTime
				and self.halfmaxGrowthNutrientConc 	== other.halfmaxGrowthNutrientConc
				and self.nutrientConsumptionRate 	== other.nutrientConsumptionRate	
				and self.marker 					== other.marker)


##################################################
##################################################


class Inoculum(object):

	def __init__(self, strain, cell_count, growth_rate=0, growth_phase='stationary', growth_cycle_timer=0):

		#------------------------------
		# Parameters:
		#------------------------------
		self.strain 			= strain

		#------------------------------
		# Variables:
		#------------------------------
		self.cellCount 			= cell_count			# N_i(t)

		self.growthRate 		= growth_rate 			# r_i(t)

		self.growthPhase 		= growth_phase 			# 'stationary', 'lag', 'exponential'
		self.growthCycleTimer	= growth_cycle_timer


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def sample(self, sample_proportion):
		try:
			sampledCellCount 	= numpy.random.poisson( self.cellCount * sample_proportion )	# num cells sampled = draw from poisson distribution given expected num cells per sample volume
		except ValueError:
			sampledCellCount 	= 0
		
		self.cellCount 			-= sampledCellCount
		
		sampledInoculum				= copy.deepcopy(self)
		sampledInoculum.cellCount 	= sampledCellCount

		return sampledInoculum


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def add(self, added_inoculum):
		assert(isinstance(added_inoculum, type(self)) and self == added_inoculum), "To add to an Inoculum, the added object must have matching type and parameters."
		self.cellCount += added_inoculum.cellCount
		return


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def __eq__(self, other):
		return(self.strain == other.strain)


##################################################
##################################################


class Plasmid(object):

	def __init__(self, name, copy_number=1, marker=''):

		self.name 				= name

		#------------------------------
		# Parameters:
		#------------------------------
		self.copyNumber 		= copy_number
		self.marker 			= marker


##################################################
##################################################


class CultureDynamics(object):

	def __init__(self):
		
		self.data 		= pandas.DataFrame(columns=['run', 'strain', 't', 't_run', 'N', 'r', 'S'])

		self.num_runs 	= 0


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def run(self, time, dt, inoculums, media, temp, noise_r=0.0, downsample_output_dt=None):

		I 	= len(inoculums)

		T 	= int(time/dt)

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# LOAD PARAMETER VALUES: 
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		rho 			= numpy.zeros(I)
		k_r 			= numpy.zeros(I)
		sigma 			= numpy.zeros(I)
		
		optimalTemp			= numpy.zeros(I)
		meanLagExitTime 	= numpy.zeros(I)
		stdevLagExitTime 	= numpy.zeros(I)

		#--------------------------------------------------

		for i, inoculum in enumerate(inoculums):

			rho[i] 			= inoculum.strain.maxGrowthRate
			k_r[i] 			= inoculum.strain.halfmaxGrowthNutrientConc
			sigma[i] 		= inoculum.strain.nutrientConsumptionRate

			optimalTemp[i] 		= inoculum.strain.optimalTemp
			meanLagExitTime[i] 	= inoculum.strain.meanLagExitTime
			stdevLagExitTime[i] = inoculum.strain.stdevLagExitTime


		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# LOAD INITIAL VARIABLE VALUES: 
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		N 		= numpy.zeros(shape=(T+1,I))	
		r 		= numpy.zeros(shape=(T+1,I))		
		S 		= numpy.zeros(shape=(T+1,1))		

		growthPhase 		= [''] * I
		growthCycleTimer	= numpy.zeros(I)

		#--------------------------------------------------

		for i, inoculum in enumerate(inoculums):

			N[0][i] 	= inoculum.cellCount
			r[0][i] 	= inoculum.growthRate

			growthPhase[i] 		= inoculum.growthPhase
			growthCycleTimer[i] = inoculum.growthCycleTimer

		S[0] = media.nutrient.concentration * media.volume


		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# RUN GROWTH MODEL:
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		for t in range(T):

			dN 			= r[t]*N[t]
			N[t+1] 		= N[t] + dN*dt
			N[t+1] 		= numpy.clip(N[t+1], a_min=0, a_max=None)

			dS 			= numpy.sum(-sigma*r[t]*N[t])
			S[t+1] 		= S[t] + dS*dt
			S[t+1] 		= numpy.clip(S[t+1], a_min=0, a_max=None)

			# Modulate the effective max growth rate by growth phase and temperature:
			r_max 		= numpy.zeros(I)
			for i in range(I):
				r_max[i] = rho[i] 					# intrinsic max exponential growth rate
				r_max[i] *= scipy.stats.norm.cdf(growthCycleTimer[i], loc=meanLagExitTime[i], scale=stdevLagExitTime[i]) 	# modeling lag phase exit (cells exit lag phase (r=0) and enter exponential phase (r=maxrate*otherfactors) at lag-exit-times that are normally distributed from cell to cell; this is approximated by the population's average growth rate following the normal distn cdf centered at the mean lag exit time
				r_max[i] *= temp/optimalTemp[i]		# growth rate dependence on temperature
				# Apply additional random noise:
				r_max[i] = numpy.random.normal(r_max[i], r_max[i]*noise_r)
			# Calculate the current growth rate as a function of nutrient availability:
			r[t+1] 		= r_max * ( S[t] / (k_r + S[t]) )
			# Update the current growth cycle phase:
			for i in range(I):
				growthCycleTimer[i] += dt
				if( (S[t] / (k_r[i] + S[t])) <= 0.01 ):
					growthPhase[i] = 'stationary'
				elif( growthPhase[i] == 'stationary' and (S[t] / (k_r[i] + S[t])) > 0.01 ):
					growthPhase[i] = 'lag'
					growthCycleTimer[i] = 0.0
				elif( growthPhase[i] == 'lag' and (rho[i]-r[t][i])/rho[i] < 0.5 ):
					growthPhase[i] = 'exponential'


			# print "N["+str(t)+"] = " + str(N[t])
			# print "r["+str(t)+"] = " + str(r[t])
			# print "S["+str(t)+"] = " + str(S[t])
			# print "------------------------------"
			# print "N["+str(t+1)+"] = " + str(N[t+1])
			# print "r["+str(t+1)+"] = " + str(r[t+1])
			# print "S["+str(t+1)+"] = " + str(S[t+1])
			# print "=============================="


		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# SAVE OUT FINAL VARIABLE VALUES: 
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		for i, inoculum in enumerate(inoculums):

			inoculum.cellCount 	= N[-1][i]
			inoculum.growthRate = r[-1][i]

		media.nutrient.concentration = S[-1] / media.volume


		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# UPDATE THE DATAFRAME HOLDING THE VARIABLE TIME SERIES:
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		t_vals 			= numpy.arange(start=0, stop=time+dt, step=dt)[:(T+1)] 			# <- last index slice is just to make sure t_vals is same length as var/param arrays generated above
		t_vals_macro 	= t_vals if self.num_runs==0 else t_vals+self.data['t'].max() 	# in terms of an overall timer, assume that each run starts immediately where the previous run left off

		for i, strainName in enumerate([inoculum.strain.name for inoculum in inoculums]):
			run_timeseries 	= {
								'run': 			[self.num_runs]*(T+1),
								'strain': 		[strainName]*(T+1),
								't': 			t_vals_macro,
								't_run':		t_vals,
								'N': 			N[:,i],
								'r': 			r[:,i],
								'S': 			S[:,0],
								'rho':			[rho[i]]*(T+1),
								'k_r':			[k_r[i]]*(T+1),
								'sigma': 		[sigma[i]]*(T+1),
								'media_volume': [media.volume]*(T+1)
							  }

			strain_run_data = pandas.DataFrame.from_dict(run_timeseries)

			self.data = self.data.append(strain_run_data)
			self.data = self.data.drop_duplicates(subset=['t', 'strain'], keep='first')


		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		self.num_runs 	+= 1

		return


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def figure(self, plot_density=True):
		import matplotlib.pyplot as pyplot
		import seaborn

		fig, ax = pyplot.subplots(2,1, figsize=(9,9))

		t_series = self.data[self.data['strain']==self.data['strain'].unique()[0]]['t'].values
		ax[0].plot(t_series, [1e6]*len(t_series), linestyle='-', color='lightgray')
		ax[1].plot(t_series, [0.0]*len(t_series), linestyle='-', color='lightgray')

		strains 		= self.data['strain'].unique()
		strain_palette 	= seaborn.color_palette("husl", len(strains))
		for i, strain in enumerate(strains):

			t_series 	= self.data[self.data['strain']==strain]['t'].values
			N_series 	= self.data[self.data['strain']==strain]['N'].values
			r_series 	= self.data[self.data['strain']==strain]['r'].values
			mediaVolume_series 	= self.data[self.data['strain']==strain]['media_volume'].values
			
			
			if(plot_density):
				ax[0].plot(t_series, N_series/mediaVolume_series, linestyle='-', color=strain_palette[i])
				ax[0].set_ylabel("population density (cells/mL)")
			else:
				ax[0].plot(t_series, N_series, linestyle='-', color=strain_palette[i])
				ax[0].set_ylabel("population size (num cells)")
			
			ax[1].plot(t_series, r_series, linestyle='-', color=strain_palette[i])
			ax[1].set_ylabel("growth rate")


		seaborn.set_style("ticks")
		seaborn.despine()

		fig.tight_layout()

		pyplot.show()


##################################################
##################################################


class Culture(object):

	def __init__(self, media=Media(volume=0), inoculums=[], dynamics=CultureDynamics(), name=""):

		self.name 		= name

		self.media		= media

		self.inoculums	= []	# class instance list has to be declared this way to avoid list being shared among instances when mutable changes are made to one instance list, apparently
		for inoculum in (inoculums if isinstance(inoculums, list) else [inoculums]):
			self.inoculums.append(inoculum)

		self.dynamics 	= dynamics


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def add(self, added_culture):
		# Add the other culture volume's media to this media:
		if(self.media is None or self.media.volume == 0.0):
			self.media = copy.deepcopy(added_culture.media)
		else:
			self.media.add(added_culture.media)
		# Add the other culture volume's inoculums to this media:
		for addedInoculum in added_culture.inoculums:
			matchingExistingInoculum = next((inoculum for inoculum in self.inoculums if inoculum == addedInoculum), None)
			if(matchingExistingInoculum):
				# A matching inoculum already exists in this culture volume, add the two inoculums together:
				matchingExistingInoculum.add(addedInoculum)
			else:
				# The added inoculum does not exist in this culture volume, add the new inoculum to the culture:
				self.inoculums.append(addedInoculum)
		# Note: Added culture volumes do not change the name or dynamics of this culture volume.
		return


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def sample(self, sample_volume, errorTolerancePct=2.0):
		originalVolume 		= self.totalVolume()
		sampledVolume 		= numpy.random.normal(sample_volume, sample_volume*(errorTolerancePct/100.0))
		sampleProportion 	= sampledVolume / originalVolume
		# Take a sample of this culture volume's media:
		sampledMedia 		= self.media.sample(sample_volume=sample_volume)	
		# Take samples of this culture volume's inoculums:
		sampledInoculums 	= []
		for inoculum in self.inoculums:
			sampledInoculums.append(inoculum.sample(sample_proportion=sampleProportion))
		# Create a new culture volume for the sample of this culture volume:
		sampledCulture 		= Culture(media=sampledMedia, inoculums=sampledInoculums, dynamics=self.dynamics, name=self.name+" sample")
		return sampledCulture


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def incubate(self, time, dt, temp):
		if(self.dynamics):
			self.dynamics.run(time=time, dt=dt, inoculums=self.inoculums, media=self.media, temp=temp)


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def totalVolume(self):
		# return numpy.sum([m.volume for m in self.media]) if len(self.media)>0 else 0.0
		return self.media.volume


	# def totalCarryingCapacity(self):
	# 	return numpy.sum([m.carryingCapacity*(m.volume/self.totalVolume()) for m in self.media]) if len(self.media)>0 else 0.0


	def getCellCounts(self):
		cellCounts = {}
		for inoculum in self.inoculums:
			cellCounts[inoculum.strain.name] = inoculum.cellCount
		return cellCounts


	def getPlasmidCounts(self):
		plasmidCounts = {}
		for inoculum in self.inoculums:
			for plasmid in inoculum.strain.plasmids:
				try:
					plasmidCounts[plasmid.name] += inoculum.cellCount * plasmid.copyNumber
				except KeyError:
					plasmidCounts[plasmid.name] = inoculum.cellCount * plasmid.copyNumber
		return plasmidCounts


	def getPlasmidFrequencies(self):
		plasmidFreqs 		= {}
		plasmidCounts 		= self.getPlasmidCounts()
		totalPlasmidCount 	= numpy.sum([plasmidCounts[plasmid] for plasmid in plasmidCounts.keys()])
		for plasmid, count in plasmidCounts.iteritems():
			plasmidFreqs[plasmid] = count/totalPlasmidCount
		return plasmidFreqs


	def getInoculumDensities(self):
		inoculumDensities = {}
		for inoculum in self.inoculums:
			inoculumDensities[inoculum.strain.name] = inoculum.cellCount / self.totalVolume()
		return inoculumDensities


	def getInoculumFrequencies(self):
		frequencies = {}
		for inoculum, count in self.getCellCounts().iteritems():
			frequencies[inoculum] = count / self.totalCellCount()
		return frequencies


	def getMarkerCounts(self):
		markerCounts = {}
		for inoculum in self.inoculums:
			markerCounts[inoculum.strain.marker] = markerCounts.get(inoculum.strain.marker, 0) + inoculum.cellCount
		return markerCounts


	def getMarkerDensities(self):
		markerDensities = {}
		markerCounts = self.getMarkerCounts()
		for marker, markerCount in markerCounts.iteritems():
			markerDensities[marker] = markerCount / self.totalVolume()
		return markerDensities


	def getSoluteConcentrations(self):
		soluteConcentrations = {}
		for solute in self.media.solutes:
			soluteConcentrations[solute.name] = solute.concentration
		return soluteConcentrations


	def getDrugConcentrations(self):
		drugConcentrations = {}
		for drug in self.media.drugs:
			drugConcentrations[drug.name] = drug.concentration
		return drugConcentrations


	def getNutrientConcentrations(self):
		nutrientConcentrations = {}
		for nutrient in self.media.nutrients:
			nutrientConcentrations[nutrient.name] = nutrient.concentration
		return nutrientConcentrations


	def totalCellCount(self):
		return numpy.sum([i.cellCount for i in self.inoculums])


	def totalCellDensity(self):
		return numpy.sum([i.cellCount/self.totalVolume() for i in self.inoculums])


	def info(self):
		print "Culture " + self.name + ":"
		print "--------------------------------"
		print "\tVolume\t\t\t= " + ("%.3f" % self.totalVolume()) + " mL"
		print "\tDensity\t\t\t= " + ( '%.2E' % Decimal(str( self.totalCellDensity() )) ) + " cfu/mL"
		print "\tTotal CellCount\t= " + ( '%.2E' % Decimal(str( self.totalCellCount() )) ) + " cfu"
		#----------
		if(len(self.inoculums)>0):
			inoculumDensities 	= self.getInoculumDensities()
			for inoculum in self.inoculums:
				print "\tInoculum " + inoculum.strain.name + ":\n" + "\t\tdensity\t\t= " + ( '%.2E' % Decimal(str( inoculumDensities[inoculum.strain.name] )) )+" cfu/mL" + "\n\t\tcell count\t= " + ( '%.2E' % Decimal(str( inoculum.cellCount )) )+" cfu" + "\n\t\tfrequency\t= " + ( '%.4f' % self.getInoculumFrequencies()[inoculum.strain.name] ) + "" + "\n\t\tgrowthPhase\t= " + inoculum.growthPhase
		else:
			print "\t(no inoculums)"
		#----------
		self.media.info()


##################################################
##################################################


class Incubator(object):

	def __init__(self, set_temp, temp_std_batch=0, temp_std_location=0, temp_std_transient=0):

		self.setTemp 			= set_temp
		self.tempSTD_batch		= temp_std_batch
		self.tempSTD_location 	= temp_std_location
		self.tempSTD_transient 	= temp_std_transient

	def incubate(self, cultures, time, dt, record_growth_curves=False):

		if(not isinstance(cultures,(list, numpy.ndarray))):
			cultures = [cultures]

		baseTemp_batch 		= numpy.random.normal(self.setTemp, self.tempSTD_batch) 

		baseTemp_locations 	= [ numpy.random.normal(baseTemp_batch, self.tempSTD_location) for culture in cultures ]

		effectiveTemps 		= baseTemp_locations

		growthCurves 	= []

		# If there is no transient temperature fluctuations 
		# (batch and location differences notwithstanding):
		if(self.tempSTD_transient == 0):
			for c, culture in enumerate(cultures):
				culture.incubate(time=time, dt=dt, temp=effectiveTemps[c])
		
		# Else, model random fluctuations in temp every dt step:
		else:
			for t in numpy.arange(start=0, stop=time, step=dt):
				
				effectiveTemps 	+= numpy.random.normal(0, self.tempSTD_transient, len(effectiveTemps))
				
				for c, culture in enumerate(cultures):
					culture.incubate(time=dt, dt=dt, temp=effectiveTemps[c])

					print "N["+str(t)+"] = " + str(culture.inoculums[0].cellCount)
					print "r["+str(t)+"] = " + str(culture.inoculums[0].growthRate)
					print "phase["+str(t)+"] = " + str(culture.inoculums[0].growthPhase)
					print "timer["+str(t)+"] = " + str(culture.inoculums[0].growthCycleTimer)
					print "END INCUBATOR T = "+str(t)
					print "------------------------------"
					if(culture.media.nutrient.concentration <= 1e-1):
						culture.info()
						exit()
					# exit()


##################################################
##################################################


class DilutionSeries(object):

	def __init__(self, initial_culture, fold_dilution=10, num_serial_dilutions=10, transfer_volume=0.030):

		self.initialCulture 	= initial_culture
		self.foldDilution 		= fold_dilution
		self.numSerialDilutions = num_serial_dilutions
		self.transferVolume 	= transfer_volume
		
		self.dilutions 	 		= [None for d in range(self.numSerialDilutions+1)]
		self.dilutions[0] 		= self.initialCulture

		for d in range(1, self.numSerialDilutions+1):
			self.dilutions[d] = Culture( 	media=[Media(volume=(self.transferVolume*self.foldDilution)-self.transferVolume, carrying_capacity=0)],
												name=self.initialCulture.name+"Serial Dilution -" + str(d)
											 )
			self.dilutions[d].add(self.dilutions[d-1].sample(self.transferVolume))

			
	def getDilution(self, dilution):
		d = numpy.abs(dilution)
		return self.dilutions[d]


	def getDensestCountableDilution(self):
		for d in range(self.numSerialDilutions):
			if(self.dilutions[d].totalCellCount()<400):
				return {'culture': self.dilutions[d], 'df': d}



##################################################
##################################################


class AgarPlate(object):

	def __init__(self, drugs=[], num_colony_landing_spots=10000, name=""):

		self.name 		= name
		self.inoculums	= []
		self.drugs 		= drugs 	# TODO: handle drug concs and effect on colonies forming baesd on cell MICs

		self.numColonyLandingSpots 	= num_colony_landing_spots
		self.colonyLandingSpots 	= [[] for n in range(num_colony_landing_spots)]


	def inoculate(self, added_culture):
		for addedInoculum in added_culture.inoculums:
			existingInoculum = next((i for i in self.inoculums if i.strain == addedInoculum.strain.name), None)
			if(existingInoculum):
				existingInoculum.cellCount += addedInoculum.cellCount
			else:
				self.inoculums.append(addedInoculum)

		for inoculum in self.inoculums:
			for cell in range(inoculum.cellCount):
				landingSpotIndex = numpy.random.randint(low=0, high=self.numColonyLandingSpots)
				if(not self.colonyLandingSpots[landingSpotIndex]):
					self.colonyLandingSpots[ landingSpotIndex ] = [inoculum.strain.marker]
				else:
					self.colonyLandingSpots[ landingSpotIndex ].append(inoculum.strain.marker)

		for l, landingSpotMarkers in enumerate(self.colonyLandingSpots):
				if(len(landingSpotMarkers) >= 1):
					self.colonyLandingSpots[l] = random.choice(landingSpotMarkers) # Of all the cells that landed on the same spot, one random one will become a visible colony
				else:
					self.colonyLandingSpots[l] = None


	def incubate(self, temp=37.0, dt=1):
		if(temp > 30):
			for l, landingSpotMarkers in enumerate(self.colonyLandingSpots):
				if(landingSpotMarkers != 0 and landingSpotMarkers != ''):
					self.colonyLandingSpots[l] = random.choice(landingSpotMarkers) # Of all the cells that landed on the same spot, one random one will become a visible colony


	def getCellCounts(self):
		cellCounts = {}
		for inoculum in self.inoculums:
			cellCounts[inoculum.strain.name] = inoculum.cellCount


	def totalCellCount(self):
		return numpy.sum([i.cellCount for i in self.inoculums])


	def getColonyCounts(self):
		colonyCounts = {}
		for marker in numpy.unique([inoculum.strain.marker for inoculum in self.inoculums]):
			colonyCounts[marker] = sum((L!=None and marker in L) for L in self.colonyLandingSpots)
		return colonyCounts


	def totalColonyCount(self):
		return sum(L!=None for L in self.colonyLandingSpots)


	def info(self):
		print "AgarPlate " + self.name + ":"
		print "--------------------------------"
		print "\tTotal Colonies = " + str(self.totalColonyCount())
		for marker in numpy.unique([inoculum.strain.marker for inoculum in self.inoculums]):
			print "\tMarker " + marker + ":\n" + "\t\tcell count\t= " + str( self.getColonyCounts()[marker] )
		# for inoculum in self.inoculums:
		# 	print "\tInoculum " + inoculum.strain.name + ":\n" + "\t\tcell count\t= " + str( self.getColonyCounts()[inoculum.strain.marker] )


##################################################
##################################################


class FlowCytometer(object):

	def __init__(self, name=""):

		self.name 		= name
		
		# self.eventSlotsPerUl 	= event_slots_per_ul
		# self.eventSlots 	= [[] for n in range(num_colony_landing_spots)]


	def read(self, sample, read_volume, event_slots_per_ul=10000):
		numEventSlots 	= event_slots_per_ul * int(read_volume*1000)
		eventSlots 		= [[] for n in range(numEventSlots)]
		for inoculum in sample.inoculums:
			for cell in range(inoculum.cellCount):
				eventSlotIndex = numpy.random.randint(low=0, high=numEventSlots)
				if(not eventSlots[eventSlotIndex]):
					eventSlots[ eventSlotIndex ] = [inoculum.strain.marker]
				else:
					eventSlots[ eventSlotIndex ].append(inoculum.strain.marker)
		cellEventCounts = {}
		for eventMarkers in eventSlots:
			overallEventMarker = '-'.join(numpy.unique(eventMarkers))
			if overallEventMarker in cellEventCounts:
				cellEventCounts[overallEventMarker] += 1
			else:
				cellEventCounts[overallEventMarker] = 1
		try:
			del cellEventCounts['']
		except KeyError:
			pass
		return cellEventCounts


		

		




















