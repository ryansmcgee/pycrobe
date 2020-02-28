from __future__ import division
from decimal import Decimal

import random
import numpy as numpy
import pandas as pandas
import scipy.stats


##################################################
##################################################


class Media:

	def __init__(self, volume, drugs=[], solutes=[], nutrient=None, isIntracellular=False):

		self.volume		= volume

		self.drugs		= []	# class instance list has to be declared this way to avoid list being shared among instances when mutable changes are made to one instance list, apparently
		for drug in (drugs if isinstance(drugs, list) else [drugs]):
			self.drugs.append(drug)

		self.solutes	= []	# class instance list has to be declared this way to avoid list being shared among instances when mutable changes are made to one instance list, apparently
		for solute in (solutes if isinstance(solutes, list) else [solutes]):
			self.solutes.append(solute)					

		self.nutrient	= nutrient

		self.isIntracellular = isIntracellular


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def sample(self, sample_proportion):
		# TODO
		return


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def add(self, added_media):
		# TODO
		return


##################################################
##################################################


class Solute:

	def __init__(self, name="Solute", concentration=0):

		self.name 			= name

		self.concentration 	= concentration


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def sample(self, sample_proportion):
		# TODO
		return


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def add(self, added_solute):
		# TODO
		return


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def __eq__(self, other):
		return(	self.name == other.name )


##################################################
##################################################


class Nutrient(Solute):

	def __init__(self, name="Nutrient", concentration=0):

		super.__init__(name, concentration)


##################################################
##################################################


class Drug(Solute): 

	def __init__(self, name="Drug", concentration=0, decay_rate=0):

		super.__init__(name, concentration)

		self.decayRate		= decay_rate									# delta_A	


##################################################
##################################################


class BetaLactam(Drug): 

	# TODO: make Drug super class, make this subclass call super

	def __init__(self, name="Beta-lactam", concentration=0, decay_rate=0, is_intracellular=False):

		super.__init__(name, concentration, decay_rate)	

		self.isIntracellular = is_intracellular				


##################################################
##################################################


class BetaLactamase(Solute):

	def __init__(self, name, decay_rate_intracellular, decay_rate_extracellular, max_hydrolysis_rate, halfmax_hydrolysis_drug_conc, is_intracellular, concentration=0):

		super.__init__(name, concentration)

		self.decayRate_intracellular 	= decay_rate_intracellular		# delta_B
		self.decayRate_extracellular 	= decay_rate_extracellular		# delta_Bext

		self.maxHydrolysisRate 	 		= max_hydrolysis_rate			# Vmax_B_i
		self.halfmaxHydrolysisDrugConc 	= halfmax_hydrolysis_drug_conc	# k_B_i

		self.isIntracellular	 		= is_intracellular				# True = in periplasm, False = in extracellular media

		self.hydrolysisRate 			= 0								# d_B_i


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def __eq__(self, other):
		return(	self.name == other.name 
				and self.decayRate_intracellular 	== other.decayRate_intracellular
				and self.decayRate_extracellular	== other.decayRate_extracellular
				and self.maxHydrolysisRate 			== other.maxHydrolysisRate
				and self.halfmaxHydrolysisDrugConc 	== other.halfmaxHydrolysisDrugConc )

		
##################################################
##################################################


class Strain:

	def __init__(self,name, max_growth_rate, optimal_temp=37.0, mean_lag_exit_time=1.5, stdev_lag_exit_time=0.125, 
							halfmax_growth_nutrient_conc=numpy.inf, nutrient_consumption_rate=0,
							marker='', plasmids=[]):

		self.name 								= name

		self.maxGrowthRate 						= max_growth_rate				# rho
		self.optimalTemp						= optimal_temp
		self.meanLagExitTime 					= mean_lag_exit_time
		self.stdevLagExitTime 					= stdev_lag_exit_time

		self.halfmaxGrowthNutrientConc 			= halfmax_growth_nutrient_conc 	# k_r		
		self.nutrientConsumptionRate 			= nutrient_consumption_rate		# sigma

		self.marker 							= marker

		self.plasmids 							= plasmids


##################################################
##################################################


class BlaStrain(Strain):

	def __init__(self,name, max_growth_rate, max_lysis_rate, halfmax_lysis_drug_conc, lysis_hill_coefficient,
							betalactamase,
							bla_production_rate, bla_leak_rate, bla_debris_sink_fraction,
							drug_diffusion_rate, drug_debris_sink_fraction,
							periplasm_volume, nutrient_consumption_rate,
							optimal_temp=37.0, mean_lag_exit_time=1.5, stdev_lag_exit_time=0.125, 
							halfmax_growth_nutrient_conc=numpy.inf, nutrient_consumption_rate=0,
							marker='', plasmids=[]):

		super().__init__(name, max_growth_rate, optimal_temp, mean_lag_exit_time, stdev_lag_exit_time, 
								halfmax_growth_nutrient_conc, nutrient_consumption_rate, marker, plasmids)

		self.maxLysisRate						= max_lysis_rate				# lamda
		self.halfmaxLysisDrugConc				= halfmax_lysis_drug_conc		# k_l
		self.lysisHillCoefficient				= lysis_hill_coefficient		# eta
		
		self.betalactamase 						= betalactamase					# BetaLactamase object

		self.betalactamaseProductionRate 		= bla_production_rate			# alpha_beta
		self.betalactamaseLeakRate 				= bla_leak_rate					# epsilon_B
		self.betalactamaseDebrisSinkFraction 	= bla_debris_sink_fraction		# xi_B

		self.betalactamDiffusionRate 	 		= drug_diffusion_rate			# epsilon_A
		self.betalactamDebrisSinkFraction  		= drug_debris_sink_fraction		# xi_A

		self.periplasmVolume 					= periplasm_volume				# Vol # in units used for media volume


##################################################
##################################################


class Inoculum:

	def __init__(self, strain, cell_count, growth_rate=0, growth_phase='stationary', growth_cycle_timer=0):

		self.strain 			= strain

		self.cellCount 			= cell_count			# N_i(t)

		self.growthRate 		= growth_rate 			# r_i(t)

		self.growthPhase 		= growth_phase 			# 'stationary', 'lag', 'exponential'
		self.growthCycleTimer	= growth_cycle_timer


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def sample(self, sample_proportion):
		# TODO
		return


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def add(self, added_inoculum):
		# TODO
		return


##################################################
##################################################


class BlaInoculum(Inoculum):

	def __init__(self, strain, cell_count, periplasm=None, growth_rate=0, lysis_rate=0, growth_phase='stationary', growth_cycle_timer=0):

		super().__init__(strain, cell_count, growth_rate, growth_phase, growth_cycle_timer)

		self.lysisRate 			= lysis_rate 			# l_i(t)

		self.periplasm 			= periplasm if periplasm is not None else Media(volume=self.strain.periplasmVolume)				# holds B_i(t) and A_i(t)


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def sample(self, sample_proportion):
		# TODO
		return


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def add(self, added_inoculum):
		# TODO
		return


##################################################
##################################################


class Plasmid:

	def __init__(self, name, copy_number=1, marker=''):

		self.name 				= name

		self.copyNumber 		= copy_number
		self.marker 			= marker


##################################################
##################################################


class CultureVolume:

	def __init__(self, media=[], inoculums=[], dynamics=CultureDynamics(), name=""):

		self.name 		= name

		self.media 		= []	# class instance list has to be declared this way to avoid list being shared among instances when mutable changes are made to one instance list, apparently
		for media in (media if isinstance(media, list) else [media]):
			self.media.append(media)

		self.inoculums	= []	# class instance list has to be declared this way to avoid list being shared among instances when mutable changes are made to one instance list, apparently
		for inoculum in (inoculums if isinstance(inoculums, list) else [inoculums]):
		# 	self.inoculums.append(inoculum)

		self.dynamics 	= dynamics


	def totalVolume(self):
		return numpy.sum([m.volume for m in self.media]) if len(self.media)>0 else 0.0


	def totalCarryingCapacity(self):
		return numpy.sum([m.carryingCapacity*(m.volume/self.totalVolume()) for m in self.media]) if len(self.media)>0 else 0.0


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


	def getCellDensities(self):
		cellDensities = {}
		for inoculum in self.inoculums:
			cellDensities[inoculum.strain.name] = inoculum.cellCount / self.totalVolume()
		return cellDensities


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
		# print markerCounts
		# exit()
		for marker, markerCount in markerCounts.iteritems():
			markerDensities[marker] = markerCount / self.totalVolume()
		return markerDensities


	def getDrugConcentrations(self):
		drugConcentrations = {}
		for drug in self.drugs:
			drugConcentrations[drug.name] = drug.mg / self.totalVolume()
		return drugConcentrations


	def totalCellCount(self):
		return numpy.sum([i.cellCount for i in self.inoculums])


	def totalCellDensity(self):
		return numpy.sum([i.cellCount/self.totalVolume() for i in self.inoculums])


	def add(self, addedCultureVolume):
		for addedMedia in addedCultureVolume.media:
			self.media.append(addedMedia)
		for addedInoculum in addedCultureVolume.inoculums:
			inoculumExistsInCulture = False
			for i, inoculum in enumerate(self.inoculums):
				if(inoculum.strain.name == addedInoculum.strain.name):
					self.inoculums[i].cellCount += addedInoculum.cellCount
					inoculumExistsInCulture = True
					break
			if(not inoculumExistsInCulture):
				self.inoculums.append(addedInoculum)
		for addedDrug in addedCultureVolume.drugs:
			drugExistsInCulture = False
			for i, drug in enumerate(self.drugs):
				if(drug.name == addedDrug.name):
					self.drugs[i].mg += addedDrug.mg
					drugExistsInCulture = True
					break
			if(not drugExistsInCulture):
				self.drugs.append(addedDrug)	


	def sample(self, volume, randomError=True, errorTolerancePct=2.0):
		presampleVolume = self.totalVolume()
		sampledVolume 	= numpy.random.normal(volume, volume*(errorTolerancePct/100.0))
		sampledMedia 	= []
		for media in self.media:
			sampledMediaVolume 	= media.volume * sampledVolume/presampleVolume
			media.volume 		-= sampledMediaVolume
			sampledMedia.append( Media(volume=sampledMediaVolume, carrying_capacity=media.carryingCapacity) )
		sampledInoculums = []
		for inoculum in self.inoculums:
			try:
				sampledCellCount 	= numpy.random.poisson( inoculum.cellCount * sampledVolume/presampleVolume )	# inoculum.cellCount * sampledVolume/presampleVolume )	# num cells sampled = draw from poisson distribution given expected num cells per sample volume
			except ValueError:
				sampledCellCount 	= 0
			inoculum.cellCount 		-= sampledCellCount
			sampledInoculums.append( Inoculum(strain=inoculum.strain, cell_count=sampledCellCount, growth_phase=inoculum.growthPhase, growth_cycle_timer=inoculum.growthCycleTimer) )
		sampledDrugs 	= []
		for drug in self.drugs:
			try:
				sampledDrugNg 	= numpy.random.poisson( drug.mg*1000000 * sampledVolume/presampleVolume )	# drug.cellCount*1000000 * sampledVolume/presampleVolume )	# num *ng* sampled = draw from poisson distribution given expected num ng's per sample volume
				sampledDrugMg	= sampledDrugNg/1000000
			except ValueError:
				sampledDrugMg 	= 0
			drug.mg 			-= sampledDrugMg
			sampledDrugs.append( Drug(name=drug.name, mg=sampledDrugMg) )
		# self.totalVolume() 		-= sampledVolume
		sampledCulture 		= CultureVolume(media=sampledMedia, inoculums=sampledInoculums, drugs=sampledDrugs, name=self.name+" sample")
		return sampledCulture


	def incubate(self, time, dt, temp):
		self.dynamics.run(time=time, dt=dt, temp=temp, inoculums=self.inoculums, media=self.media)


	def info(self):
		print "CultureVolume " + self.name + ":"
		print "--------------------------------"
		print "\tVolume\t\t\t= " + ("%.3f" % self.totalVolume()) + " mL"
		print "\tCarrying Cap.\t= " + ( '%.2E' % Decimal(str( self.totalCarryingCapacity() )) ) + " cfu/mL"
		print "\tDensity\t\t\t= " + ( '%.2E' % Decimal(str( self.totalCellDensity() )) ) + " cfu/mL"
		print "\tTotal CellCount\t= " + ( '%.2E' % Decimal(str( self.totalCellCount() )) ) + " cfu"
		if(len(self.inoculums)>0):
			for inoculum in self.inoculums:
				print "\tInoculum " + inoculum.strain.name + ":\n" + "\t\tdensity\t\t= " + ( '%.2E' % Decimal(str( inoculum.cellCount/self.totalVolume() )) )+" cfu/mL" + "\n\t\tcell count\t= " + ( '%.2E' % Decimal(str( inoculum.cellCount )) )+" cfu" + "\n\t\tgrowthPhase\t= " + inoculum.growthPhase
		else:
			print "\t(no inoculums)"
		if(len(self.drugs)>0):
			for drug in self.drugs:
				print "\tDrug " + drug.name + ":\n" + "\t\tconcentration\t= " + ("%.5f" % self.getDrugConcentrations()[drug.name]) + " mg/mL"
		else:
				print "\t(no drugs)"


##################################################
##################################################


class DilutionSeries:

	def __init__(self, initial_culture, fold_dilution=10, num_serial_dilutions=10, transfer_volume=0.030):

		self.initialCulture 	= initial_culture
		self.foldDilution 		= fold_dilution
		self.numSerialDilutions = num_serial_dilutions
		self.transferVolume 	= transfer_volume
		
		self.dilutions 	 		= [None for d in range(self.numSerialDilutions+1)]
		self.dilutions[0] 		= self.initialCulture

		for d in range(1, self.numSerialDilutions+1):
			self.dilutions[d] = CultureVolume( 	media=[Media(volume=(self.transferVolume*self.foldDilution)-self.transferVolume, carrying_capacity=0)],
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


class AgarPlate:

	def __init__(self, drugs=[], num_colony_landing_spots=10000, name=""):

		self.name 		= name
		self.inoculums	= []
		self.drugs 		= drugs 	# TODO: handle drug concs and effect on colonies forming baesd on cell MICs

		self.numColonyLandingSpots 	= num_colony_landing_spots
		self.colonyLandingSpots 	= [[] for n in range(num_colony_landing_spots)]


	def inoculate(self, addedCultureVolume):
		for addedInoculum in addedCultureVolume.inoculums:
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


	def incubate(self, dt=1, temp=37.0):
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


class Incubator:

	def __init__(self, set_temp, temp_std_batch=0, temp_std_location=0, temp_std_transient=0):

		self.setTemp 			= set_temp
		self.tempSTD_batch		= temp_std_batch
		self.tempSTD_location 	= temp_std_location
		self.tempSTD_transient 	= temp_std_transient

	def incubate(self, dt, cultures, time, record_growth_curves=False):

		if(not isinstance(cultures,(list, numpy.ndarray))):
			cultures = [cultures]

		baseTemp_batch 		= numpy.random.normal(self.setTemp, self.tempSTD_batch) 

		baseTemp_locations 	= [ numpy.random.normal(baseTemp_batch, self.tempSTD_location) for culture in cultures ]

		effectiveTemps 		= baseTemp_locations

		growthCurves 	= []

		t=0
		for t in numpy.arange(start=0, stop=time, step=dt):
			
			effectiveTemps 	+= numpy.random.normal(0, self.tempSTD_transient, len(effectiveTemps))
			
			for c, culture in enumerate(cultures):
				culture.incubate(dt=dt, temp=effectiveTemps[c])

				if(record_growth_curves):
					growthCurves.append( {'culture':culture.name, 'time':t+dt, 'temp':effectiveTemps[c], 'density':culture.totalCellDensity()} )

		if(record_growth_curves):
			return pandas.DataFrame(growthCurves)




##################################################
##################################################


class FlowCytometer:

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


		

		




















