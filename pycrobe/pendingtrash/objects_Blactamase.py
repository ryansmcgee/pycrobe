from __future__ import division
from decimal import Decimal

import random
import numpy as numpy
import pandas as pandas
import scipy.stats


##################################################
##################################################


class Media:#_BetaLactamase:

	def __init__(self, volume, drugs=[], solutes=[], nutrient=None, isIntracellular=False):

		self.volume				= volume

		self.drugs	= []	# class instance list has to be declared this way to avoid list being shared among instances when mutable changes are made to one instance list, apparently
		for drug in (drugs if isinstance(drugs, list) else [drugs]):
			self.drugs.append(drug)

		self.solutes	= []	# class instance list has to be declared this way to avoid list being shared among instances when mutable changes are made to one instance list, apparently
		for solute in (solutes if isinstance(solutes, list) else [solutes]):
			self.solutes.append(solute)					

		self.nutrient 			= nutrient

		self.isIntracellular	= isIntracellular

	
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
			
	def incubate(self, dt, temp, inoculums):
		
		for drug in self.drugs:
			drug.update(dt, temp, inoculums, self.isIntracellular)

		for solute in self.solutes:
			solute.update(dt, temp, inoculums, self.isIntracellular)

		self.nutrient.update(dt, temp, inoculums, self.isIntracellular)


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def sample(self, sample_proportion):
		return


##################################################
##################################################


class Drug_BetaLactam:

	def __init__(self, name, concentration=0, decayRate=0):

		self.name 			= name

		self.concentration 	= concentration

		self.decayRate		= decayRate					^ delta_A	


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def update(self, dt, temp, inoculums, isIntracellular=False):

		if(self.isIntracellular): 						^ A_in(t)

			# Assume the inoculum for which this beta-lactam drug
			# is intracellular is the only inoculum in the list:
			inoculum_i 	= inoculums[0]

			A_i 		= self.concentration
			delta_A 	= self.decayRate
			
			# Note: Diffusion term is handled in Inoculum_BetaLactamase.incubate()

			d_beta_i 	= self.strain.

			# Find the beta-lactamase solute in the inoculum for which this B-lactam drug is intracellular:
			betalactamase_i = None
			for solute in inoculum_i.periplasmMedia.solutes:
				if(isinstance(solute, BetaLactamase)):
					betalactamase_i = solute

			if(betalactamase_i is not None):
				beta_i 		= betalactamase_i.concentration
				d_beta_i	= betalactamase_i.degradationRate_intracellular

				dA_i = -delta_A*A_i - d_B_i*beta_i 

			else:
				# No intracellular beta-lactamase found in this inoculum.
				dA_i 	= -delta_A*A_i 

			#------------------------------

			self.concentration += dA_i*dt

		#------------------------------

		else: #isExtracellular							^ A_ext(t)

			A_ext 		= self.concentration
			delta_A 	= self.decayRate

			degradationTerm = 0
			for solute in 
			












			for inoculum_i in inoculums:
				B_i 	= inoculum_i
				d_B_i 	= inoculum_i.degradationRate_extracellular

				# Find the beta-lactamase solute in this inoculum_i:
				betalactamase_i = None
				for solute in inoculum_i.periplasmMedia.solutes:
					if(isinstance(solute, BetaLactamase)):
						betalactamase_i = solute

				if(betalactamase_i is not None):
					beta_i 		= betalactamase_i.concentration
					d_beta_i	= betalactamase_i.degradationRate_intracellular

					dA_i = -delta_A*A_i - d_B_i*beta_i 

				else:
					# No intracellular beta-lactamase found in this inoculum.
					dA_i 	= -delta_A*A_i 

				#------------------------------








		return self.concentration


##################################################
##################################################


class Nutrient:

	def __init__(self):

		self.name 				= name

		self.concentration


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def update(self, dt, temp, inoculums, isIntracellular=False):

		if(self.isIntracellular): 						^ S(t)

		else:											# not handling intracellular nutrients in this 


##################################################
##################################################


class BetaLactamase:

	def __init__(self, name,   ):

		self.name 				= name

		self.concentration

		

		# self.halfmaxLysisEnzymeConc						^ k_beta_i
		
		# self.degradationRate_intracellular 				^ d_beta
		# self.degradationRate_extracellular 				^ d_B

		# self.maxDegradationRate_intracellular 			^ Vmax_beta_i
		# self.halfmaxDegradationRateDrugConc_intracellular 	^k_beta_i

		# self.maxDegradationRate_extracellular 			^ Vmax_B_i
		# self.halfmaxDegradationRateDrugConc_extracellular 	^k_beta_i

		self.maxDegradationRate 						^ Vmax_beta_i
		self.halfmaxDegradationDrugConc 				^ k_beta_i

		self.decayRate_intracellular 					^ delta_beta
		self.decayRate_extracellular 					^ delta_B

		self.isIntracellular							# True = in periplams, False = in extracellular media


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def degradationRate(self):

		if(self.isIntracellular): 						^ d_beta_i(t)

		else:											^ d_B_i(t)

		return


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

^_^	def update(self, dt, temp, inoculums, isIntracellular=False):

		if(isIntracellular): 							^ beta_i(t)
			
			# Assume the inoculum for which this beta-lactamase 
			# is intracellular is the only inoculum in the list:
			inoculum_i 	= inoculums[0] 	

			alpha_beta 	= self.betalactamaseProductionRate
			delta_beta 	= self.decayRate_intracellular
			r_i 		= inoculums_i.growthRate
			
			beta_i 		= alpha_beta / (r_i + delta_beta)

			self.concentration 	= beta_i

		#------------------------------

		else: # isExtracellular							^ B(t)

			# Find the inoculum whose intracellular beta-lactamase is the 
			# same type as this extracellular beta-lactamase:
			inoculum_i 		= None
			betalactamase_i = None
			for inoculum in inoculums:
				for solute in inoculum.periplasmMedia.solutes:
					if(isinstance(solute, BetaLactamase) and solute.name == self.name):
						inoculum_i 		= inoculum
						betalactamase_i = solute

			B_i 		= self.concentration
			delta_B 	= self.decayRate_extracellular
			
			if(inoculum_i is not None):
				l_i 		= inoculum_i.lysisRate
				beta_i 		= betalactamase_i.concentration
				V_i			= inoculum_i.strain.periplasmVolume
				N_i 		= inoculum_i.cellCount

				dB_i 	= -delta_B*B_i + l_i*( beta_i*V_i*N_i )
			
			else:
				# No intracellular beta-lactamases the same type as 
				# this extracellular beta-lactamase could be found.
				dB_i 	= -delta_B*B_i 

			#------------------------------

			self.concentration += dB_i*dt

		#------------------------------

		return self.concentration


##################################################
##################################################


class Strain_BetaLactamase:

	def __init__(self, name, max_growth_rate, optimal_temp=37.0, mean_lag_exit_time=1.5, stdev_lag_exit_time=0.125, marker='', plasmids=[]):

		# self.name 				= name

		self.maxGrowthRate 		= max_growth_rate		^ rho
		self.halfmaxGrowthNutrientConc 				^ k_R
		
		self.optimalTemp		= optimal_temp
		self.meanLagExitTime 	= mean_lag_exit_time
		self.stdevLagExitTime 	= stdev_lag_exit_time

		self.maxLysisRate								^ lamda
		self.halfmaxLysisDrugConc						^ k_L
		self.lysisHillCoefficient						^ h
		
		self.betalactamase 								# BetaLactamase object

		self.betalactamaseProductionRate				^ alpha_beta

		self.periplasmVolume 							^ V # in units used for media volume

		self.nutrientConsumptionRate 					^ sigma

		self.membraneDiffusionRate 						^ epsilon


		# self.marker 			= marker

		# self.plasmids 			= plasmids


##################################################
##################################################


class Inoculum_BetaLactamase:

	def __init__(self, strain, cell_count, growth_phase='stationary', growth_cycle_timer=0):

		self.strain 			= strain

		self.cellCount 			= cell_count			^ N_i(t)

		self.growthRate 		= 

		self.lysisRate 			= 

		self.growthPhase 		= growth_phase # 'stationary', 'lag', 'exponential'
		self.growthCycleTimer	= growth_cycle_timer

		self.periplasmMedia 	= Media() 				# holds beta_i(t) and A_i(t)


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

^_^	def updateGrowthRate(self, dt, temp, incubationMedia, noise=0.05):	^ r_i(t)

		# Modulate the effective max growth rate by growth phase:
		r_max 	= self.strain.maxGrowthRate 		# intrinsic max exponential growth rate
		r_max 	*= scipy.stats.norm.cdf(self.growthCycleTimer, loc=self.strain.meanLagExitTime, scale=self.strain.stdevLagExitTime) 	# modeling lag phase exit (cells exit lag phase (r=0) and enter exponential phase (r=maxrate*otherfactors) at lag-exit-times that are normally distributed from cell to cell; this is approximated by the population's average growth rate following the normal distn cdf centered at the mean lag exit time

		# Modulate the effective max growth rate by incubation temperature:
		r_max 	*= temp/self.strain.optimalTemp		# growth rate dependence on temperature
		
		# Calculate the current growth rate as a function of nutrient availability:
		S 		= incubationMedia.nutrient.concentration
		k_R 	= self.strain.halfmaxGrowthNutrientConc
		r 		= r_max * ( S / (k_R + S) )

		# Apply additional random noise:
		r 		= numpy.random.normal(r, r*noise)	

		# Update the current growth cycle phase:
		self.growthCycleTimer += dt
		if(saturationCoeff <= 0.01):
			self.growthPhase = 'stationary'
		elif(self.growthPhase == 'stationary' and saturationCoeff > 0.01):
			self.growthPhase = 'lag'
			self.growthCycleTimer = 0.0
		elif(self.growthPhase == 'lag' and (self.strain.maxGrowthRate-r)/self.strain.maxGrowthRate < 0.01):
			self.growthPhase = 'exponential'

		# print str(self.strain.name) + ' ' + str(self.strain.maxGrowthRate) + ' ' + str(r)+ ' ' + str(self.growthPhase) + ' ' + str(self.growthCycleTimer)
		
		#Update growth rate variable:
		self.growthRate = r

		return self.growthRate


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

^_^	def updateLysisRate(self):							^ l_i(t)

		l_max	= self.strain.maxLysisRate
		k_L 	= self.strain.halfmaxLysisDrugConc
		A_i 	= next((x for x in self.periplasmMedia.drugs if isinstance(x, Drug_BetaLactam)), None).concentration
		h 		= self.strain.lysisHillCoefficient
		r 		= self.growthRate

		l 		= l_max * ( A_i**h / (k_L**h + A_i**h) ) * r

		# Update lysis rate variable:
		self.lysisRate = l

		return self.lysisRate


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def incubate(self, dt, temp, media):

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Update growth and lysis rates:
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		self.updateGrowthRate(dt, temp, media)
		self.updateLysisRate()

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Update cell count:
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		dN 	= ( self.growthRate - self.lysisRate ) * self.cellCount
		self.cellCount += dN*dt

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Update periplasm:
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Handle diffusion of drug in/out of the periplasm:
		drug_ext = None
		A_ext	 = 0
		for drug_e in media.drugs:	
			if(type(drug_e) == Drug_BetaLactam):
				drug_ext 	= drug_e
				A_ext 		= drug_ext.concentration
				# There should only be one drug object of each type (eg Drug_BetaLactam) in the list of drugs in a Media object (like drugs should be combined in CultureVolume add() calls.)
				break
		for drug_i in self.periplasmMedia.drugs:	
			if(type(drug_i) == Drug_BetaLactam):
				A_i 	= drug_i.concentration
				epsilon = self.strain.membraneDiffusionRate
				#------------------------------
				dA_i 	= epsilon * ( A_ext - A_i )
				drug_i.concentration += dA_i*dt

				if(drug_ext is not None):
					drug_ext.concentration -= dA_i*dt # the drug conc that enters the cell leaves the external media
				#------------------------------
				# There should only be one drug object of each type (eg Drug_BetaLactam) in the list of drugs in a Media object (like drugs should be combined in CultureVolume add() calls.)
				break
		#------------------------------
		# Update periplasm media contents:
		self.periplasmMedia.incubate(dt, temp, inoculums=[self])
	

	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def sample(self, sample_proportion):
		return


##################################################
##################################################


class Plasmid:

# 	def __init__(self, name, copy_number=1, marker=''):

# 		self.name 				= name

# 		self.copyNumber 		= copy_number
# 		self.marker 			= marker


##################################################
##################################################


class CultureVolume:

	def __init__(self, media=[], inoculums=[], drugs=[], name=""):

		# self.name 		= name

		# self.media	= []	# class instance list has to be declared this way to avoid list being shared among instances when mutable changes are made to one instance list, apparently
		# for media in (media if isinstance(media, list) else [media]):
		# 	self.media.append(media)

		# self.inoculums	= []	# class instance list has to be declared this way to avoid list being shared among instances when mutable changes are made to one instance list, apparently
		# for inoculum in (inoculums if isinstance(inoculums, list) else [inoculums]):
		# 	self.inoculums.append(inoculum)

		# self.drugs	= []	# class instance list has to be declared this way to avoid list being shared among instances when mutable changes are made to one instance list, apparently
		# for drug in (drugs if isinstance(drugs, list) else [drugs]):
		# 	self.drugs.append(drug)


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


	def incubate(self, dt, temp):
		saturationCoeff = (self.totalCarryingCapacity() - self.totalCellDensity())/self.totalCarryingCapacity()
		for inoculum in self.inoculums:
			inoculum.incubate(dt=dt, temp=temp saturationCoeff=saturationCoeff)


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


		

		




















