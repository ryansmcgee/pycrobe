from __future__ import division
from decimal import Decimal

import numpy
import pandas
import scipy.stats
import copy

from standard import *


##################################################
##################################################


class BetaLactam(Drug): 

	def __init__(self, name="Beta-lactam", concentration=0, decay_rate=0.00385, is_intracellular=False):

		super(BetaLactam, self).__init__(name, concentration, decay_rate)	

		#------------------------------
		# Parameters:
		#------------------------------
		# -> # Exponential decay parameter set to 0.00385 such that t_90 (time where conc = 90% orig conc) = 27 hrs
		#		as reported for cefatriaxone (a cephalosporin related to CTX) in diluted human serum in Esteban et al. 1990.
		self.isIntracellular = is_intracellular	


##################################################
##################################################


class BetaLactamase(Solute):

	def __init__(self, name, decay_rate_intracellular, decay_rate_extracellular, max_hydrolysis_rate, halfmax_hydrolysis_drug_conc, is_intracellular, concentration=0, hydrolysis_rate=0):

		super(BetaLactamase, self).__init__(name, concentration)

		#------------------------------
		# Parameters:
		#------------------------------
		self.decayRate_intracellular 	= decay_rate_intracellular		# delta_B
		self.decayRate_extracellular 	= decay_rate_extracellular		# delta_Bext

		self.maxHydrolysisRate 	 		= max_hydrolysis_rate			# Vmax_B_i
		self.halfmaxHydrolysisDrugConc 	= halfmax_hydrolysis_drug_conc	# k_B_i

		self.isIntracellular	 		= is_intracellular				# True = in periplasm, False = in extracellular media

		#------------------------------
		# Variables:
		#------------------------------
		self.hydrolysisRate 			= hydrolysis_rate				# d_B_i


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def __eq__(self, other):
		return(self.name == other.name 
				and type(self) 						== type(other)
				and self.decayRate_intracellular 	== other.decayRate_intracellular
				and self.decayRate_extracellular	== other.decayRate_extracellular
				and self.maxHydrolysisRate 			== other.maxHydrolysisRate
				and self.halfmaxHydrolysisDrugConc 	== other.halfmaxHydrolysisDrugConc )

		
##################################################
##################################################


class BlaStrain(Strain):

	def __init__(self,name, max_growth_rate, max_lysis_rate, halfmax_lysis_drug_conc, lysis_hill_coefficient,
							betalactamase,
							bla_production_rate, bla_saturation_conc, halfmax_bla_production_conc,
							bla_leak_rate, bla_debris_sink_fraction,
							drug_diffusion_rate, drug_debris_sink_fraction,
							periplasm_volume,
							optimal_temp=37.0, mean_lag_exit_time=1.5, stdev_lag_exit_time=0.125, 
							halfmax_growth_nutrient_conc=numpy.inf, nutrient_consumption_rate=0,
							marker='', plasmids=[]):

		super(BlaStrain, self).__init__(name, max_growth_rate, optimal_temp, mean_lag_exit_time, stdev_lag_exit_time, 
										halfmax_growth_nutrient_conc, nutrient_consumption_rate, marker, plasmids)

		#------------------------------
		# Parameters:
		#------------------------------
		self.maxLysisRate						= max_lysis_rate				# lamda
		self.halfmaxLysisDrugConc				= halfmax_lysis_drug_conc		# k_l
		self.lysisHillCoefficient				= lysis_hill_coefficient		# eta

		self.betalactamaseProductionRate 		= bla_production_rate			# alpha_beta
		self.betalactamaseSaturationConc 		= bla_saturation_conc 			# B_max
		self.halfmaxBlaProductionBlaConc 		= halfmax_bla_production_conc	# k_alpha
		self.betalactamaseLeakRate 				= bla_leak_rate					# epsilon_B
		self.betalactamaseDebrisSinkFraction 	= bla_debris_sink_fraction		# xi_B

		self.betalactamDiffusionRate 	 		= drug_diffusion_rate			# epsilon_A
		self.betalactamDebrisSinkFraction  		= drug_debris_sink_fraction		# xi_A

		self.periplasmVolume 					= periplasm_volume				# Vol # in units used for media volume

		#------------------------------
		# Variable-holding Objects:
		#------------------------------
		self.betalactamase 						= betalactamase					# BetaLactamase object

		
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
				and self.marker 					== other.marker
				and self.maxLysisRate						== other.maxLysisRate				
				and self.halfmaxLysisDrugConc				== other.halfmaxLysisDrugConc		
				and self.lysisHillCoefficient				== other.lysisHillCoefficient		
				and self.betalactamase 						== other.betalactamase					
				and self.betalactamaseProductionRate 		== other.betalactamaseProductionRate			
				and self.betalactamaseLeakRate 				== other.betalactamaseLeakRate					
				and self.betalactamaseDebrisSinkFraction 	== other.betalactamaseDebrisSinkFraction		
				and self.betalactamDiffusionRate 	 		== other.betalactamDiffusionRate			
				and self.betalactamDebrisSinkFraction  		== other.betalactamDebrisSinkFraction		
				and self.periplasmVolume 					== other.periplasmVolume
				)


##################################################
##################################################


class BlaInoculum(Inoculum):

	def __init__(self, strain, cell_count, periplasm=None, growth_rate=0, lysis_rate=0, growth_phase='stationary', growth_cycle_timer=0):

		super(BlaInoculum, self).__init__(strain, cell_count, growth_rate, growth_phase, growth_cycle_timer)

		#------------------------------
		# Variables:
		#------------------------------
		self.lysisRate 			= lysis_rate 			# l_i(t)

		#------------------------------
		# Variable-holding Objects:
		#------------------------------
		self.periplasm 			= periplasm if periplasm is not None else Media(volume=self.strain.periplasmVolume, solutes=[self.strain.betalactamase], is_intracellular=True)				# holds B_i(t) and A_i(t)


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
		self.periplasm.add(added_inoculum.periplasm)
		return


###################################################
###################################################


class BetaLactamaseDynamics(CultureDynamics):

	def __init__(self):
		super(BetaLactamaseDynamics, self).__init__()
		return


	def run(self, time, dt, inoculums, media, temp, noise_r=0.0, downsample_output_dt=None):

		I 	= len(inoculums)

		T 	= int(time/dt)

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# LOAD PARAMETER VALUES: 
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		rho 			= numpy.zeros(I)
		k_r 			= numpy.zeros(I)
		lamda 			= numpy.zeros(I)
		k_l 			= numpy.zeros(I)
		eta 			= numpy.zeros(I)
		alpha_B 		= numpy.zeros(I)
		B_max 			= numpy.zeros(I)
		k_alpha 		= numpy.zeros(I)
		delta_B 		= numpy.zeros(I)
		delta_Bext 		= numpy.zeros(I)
		epsilon_B 		= numpy.zeros(I)
		xi_B 			= numpy.zeros(I)
		delta_A 		= numpy.zeros(I)
		epsilon_A 		= numpy.zeros(I)
		xi_A 			= numpy.zeros(I)
		Vmax_B 			= numpy.zeros(I)
		k_B 			= numpy.zeros(I)
		Vmax_Bext 		= numpy.zeros(I)
		k_Bext 			= numpy.zeros(I)
		Vol 			= numpy.zeros(I)
		sigma 			= numpy.zeros(I)

		optimalTemp			= numpy.zeros(I)
		meanLagExitTime 	= numpy.zeros(I)
		stdevLagExitTime 	= numpy.zeros(I)

		#--------------------------------------------------

		for i, inoculum in enumerate(inoculums):

			rho[i] 			= inoculum.strain.maxGrowthRate
			k_r[i] 			= inoculum.strain.halfmaxGrowthNutrientConc
			lamda[i] 		= inoculum.strain.maxLysisRate
			k_l[i] 	 		= inoculum.strain.halfmaxLysisDrugConc
			eta[i] 			= inoculum.strain.lysisHillCoefficient
			alpha_B[i] 		= inoculum.strain.betalactamaseProductionRate
			B_max[i] 		= inoculum.strain.betalactamaseSaturationConc
			k_alpha[i] 		= inoculum.strain.halfmaxBlaProductionBlaConc
			delta_B[i] 		= inoculum.strain.betalactamase.decayRate_intracellular
			delta_Bext[i]	= inoculum.strain.betalactamase.decayRate_extracellular
			epsilon_B[i] 	= inoculum.strain.betalactamaseLeakRate
			xi_B[i] 		= inoculum.strain.betalactamaseDebrisSinkFraction
			epsilon_A[i] 	= inoculum.strain.betalactamDiffusionRate
			xi_A[i] 		= inoculum.strain.betalactamDebrisSinkFraction
			Vmax_B[i] 		= inoculum.strain.betalactamase.maxHydrolysisRate
			k_B[i] 			= inoculum.strain.betalactamase.halfmaxHydrolysisDrugConc
			Vmax_Bext[i] 	= inoculum.strain.betalactamase.maxHydrolysisRate
			k_Bext[i] 		= inoculum.strain.betalactamase.halfmaxHydrolysisDrugConc
			Vol[i] 			= inoculum.strain.periplasmVolume
			sigma[i] 		= inoculum.strain.nutrientConsumptionRate

			optimalTemp[i] 		= inoculum.strain.optimalTemp
			meanLagExitTime[i] 	= inoculum.strain.meanLagExitTime
			stdevLagExitTime[i] = inoculum.strain.stdevLagExitTime

		#--------------------------------------------------

		for drug in media.drugs:
			if(isinstance(drug, BetaLactam)):
				delta_A = drug.decayRate
				# Assume Media objects hold only one drug object per type (BetaLactam)
				break


		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# LOAD INITIAL VARIABLE VALUES: 
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		N 		= numpy.zeros(shape=(T+1,I))	
		r 		= numpy.zeros(shape=(T+1,I))	
		l 		= numpy.zeros(shape=(T+1,I))	
		B 		= numpy.zeros(shape=(T+1,I))	
		Bext 	= numpy.zeros(shape=(T+1,I))	
		d_B 	= numpy.zeros(shape=(T+1,I))	
		d_Bext 	= numpy.zeros(shape=(T+1,I))	
		A 		= numpy.zeros(shape=(T+1,I))	
		Aext 	= numpy.zeros(shape=(T+1,1)) 	
		S 		= numpy.zeros(shape=(T+1,1))		

		growthPhase 		= [''] * I
		growthCycleTimer	= numpy.zeros(I)

		#--------------------------------------------------

		for i, inoculum in enumerate(inoculums):

			N[0][i] 	= inoculum.cellCount
			r[0][i] 	= inoculum.growthRate
			l[0][i] 	= inoculum.lysisRate

			growthPhase[i] 		= inoculum.growthPhase
			growthCycleTimer[i] = inoculum.growthCycleTimer

			for solute in inoculum.periplasm.solutes:
				if(isinstance(solute, BetaLactamase)):
					B[0][i] 	= solute.concentration
					d_B[0][i] 	= solute.hydrolysisRate

					# Find the corresponding extracellular BetaLactamase:
					matchingBlaFound = False
					for solute_ext in media.solutes:
						if(isinstance(solute_ext, BetaLactamase) and solute == solute_ext):
							Bext[0][i] 	= solute_ext.concentration
							d_Bext[0][i] 	= solute_ext.hydrolysisRate
							matchingBlaFound = True
							break
					if(not matchingBlaFound): 
						# Make new object for Bla type that doesn't currently exist in extracellular media:
						# betalactamase_ext = BetaLactamase(name=solute.name)
						# betalactamase_ext.__dict__.update(solute.__dict__)
						betalactamase_ext 	= copy.deepcopy(solute)
						betalactamase_ext.concentration = 0.0
						betalactamase_ext.isIntracellular = False
						media.solutes.append(betalactamase_ext)

					# Assume Media objects hold only one solute object per type (BetaLactamase)
					break

			A[0][i] = 0.0
			for drug in inoculum.periplasm.drugs:
				if(isinstance(drug, BetaLactam)):
					A[0][i] 	= drug.concentration
					# Assume Media objects hold only one drug object per type (BetaLactam)
					break

		#--------------------------------------------------

		Aext[0] = 0.0
		for drug in media.drugs:
			if(isinstance(drug, BetaLactam)):
				Aext[0] = drug.concentration
				# Assume Media objects hold only one drug object per type (BetaLactam)
				break

		S[0] = media.nutrient.concentration * media.volume


		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# RUN GROWTH MODEL:
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		for t in range(T):

			dN 			= (r[t] - l[t])*N[t]
			N[t+1] 		= N[t] + dN*dt
			N[t+1] 		= numpy.clip(N[t+1], a_min=0, a_max=None)

			dB 			= -delta_B*B[t] - r[t]*B[t] - epsilon_B*(B[t] - Bext[t]) + alpha_B*( (B_max - B[t]) / (k_alpha + (B_max - B[t])) )
			B[t+1] 		= B[t] + dB*dt
			B[t+1] 		= numpy.clip(B[t+1], a_min=0, a_max=None)

			dBext 		= -delta_Bext*Bext[t] + epsilon_B*(B[t] - Bext[t])*Vol*N[t] + l[t]*(1-xi_B)*B[t]*Vol*N[t]
			Bext[t+1] 	= Bext[t] + dBext*dt
			Bext[t+1] 		= numpy.clip(Bext[t+1], a_min=0, a_max=None)

			dA 			= -delta_A*A[t] - d_B[t]*B[t] - r[t]*A[t] + epsilon_A*(Aext[t] - A[t])
			A[t+1] 		= A[t] + dA*dt
			A[t+1] 		= numpy.clip(A[t+1], a_min=0, a_max=None)

			dAext 		= -delta_A*Aext[t] - numpy.sum(d_Bext*Bext[t]) - numpy.sum(epsilon_A*(Aext[t] - A[t])*Vol*N[t]) + numpy.sum(l[t]*(1-xi_A)*A[t]*Vol*N[t])
			Aext[t+1] 	= Aext[t] + dAext*dt
			Aext[t+1] 	= numpy.clip(Aext[t+1], a_min=0, a_max=None)

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

			l[t+1] 		= 2 * lamda * ( A[t+1]**eta / (k_l**eta + A[t+1]**eta) ) * r[t+1]

			d_B[t+1]	= (Vmax_B * A[t+1])/(k_B + A[t+1])

			d_Bext[t+1]	= (Vmax_Bext * Aext[t+1])/(k_Bext + Aext[t+1])


			# print "N["+str(t)+"] = " + str(N[t])
			# print "r["+str(t)+"] = " + str(r[t])
			# print "l["+str(t)+"] = " + str(l[t])
			# print "r-l["+str(t)+"] = " + str(r[t]-l[t])
			# print "B["+str(t)+"] = " + str(B[t])
			# print str(l[t]*(1-xi_B)*B[t]*Vol*N[t])
			# print "Bext["+str(t)+"] = " + str(Bext[t])
			# print "A["+str(t)+"] = " + str(A[t])
			# print "Aext["+str(t)+"] = " + str(Aext[t])
			# print "d_B["+str(t)+"] = " + str(d_B[t])
			# print "d_Bext["+str(t)+"] = " + str(d_Bext[t])
			# print "S["+str(t)+"] = " + str(S[t])
			# print "------------------------------"
			# print "N["+str(t+1)+"] = " + str(N[t+1])
			# print "r["+str(t+1)+"] = " + str(r[t+1])
			# print "l["+str(t+1)+"] = " + str(l[t+1])
			# print "B["+str(t+1)+"] = " + str(B[t+1])
			# print "Bext["+str(t+1)+"] = " + str(Bext[t+1])
			# print "A["+str(t+1)+"] = " + str(A[t+1])
			# print "Aext["+str(t+1)+"] = " + str(Aext[t+1])
			# print "d_B["+str(t+1)+"] = " + str(d_B[t+1])
			# print "d_Bext["+str(t+1)+"] = " + str(d_Bext[t+1])
			# print "S["+str(t+1)+"] = " + str(S[t+1])
			# print "=============================="


		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# SAVE OUT FINAL VARIABLE VALUES: 
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		for i, inoculum in enumerate(inoculums):

			inoculum.cellCount 	= N[-1][i]
			inoculum.growthRate = r[-1][i]
			inoculum.lysisRate 	= l[-1][i]

			for solute in inoculum.periplasm.solutes:
				if(isinstance(solute, BetaLactamase)):
					solute.concentration 	= B[-1][i]
					solute.hydrolysisRate 	= d_B[-1][i]

					# Update the corresponding extracellular BetaLactamase:
					matchingExtBla = next((solute_ext for solute_ext in media.solutes if solute == solute_ext), None)
					if(matchingExtBla):
						# A matching extracellular betalactamase already exists in this culture, 
						# update the extracellular betalactamase variables:
						matchingExtBla.concentration 	= Bext[-1][i]
						matchingExtBla.hydrolysisRate 	= d_Bext[-1][i]
					else:
						# The added extracellular betalactamase does not exist in this culture, 
						# add the new extracellular betalactamase to the culture:
						newExtBla = copy.deepcopy(solute)
						newExtBla.concentration 	= Bext[-1][i]
						newExtBla.hydrolysisRate 	= d_Bext[-1][i]
						media.solutes.append(newExtBla)
					# Assume Media objects hold only one solute object per type (BetaLactamase)
					break

			for drug in inoculum.periplasm.drugs:
				if(isinstance(drug, BetaLactam)):
					drug.concentration 	= A[-1][i]
					# Assume Media objects hold only one drug object per type (BetaLactam)
					break

		#--------------------------------------------------

		for drug in media.drugs:
			if(isinstance(drug, BetaLactam)):
				drug.concentration = Aext[-1]
				# Assume Media objects hold only one drug object per type (BetaLactam)
				break

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
								'l': 			l[:,i],
								'B': 			B[:,i],
								'Bext': 		Bext[:,i],
								'd_B': 			d_B[:,i],
								'd_Bext': 		d_Bext[:,i],
								'A': 			A[:,i],
								'Aext': 		Aext[:,0],
								'S': 			S[:,0],
								'rho':			[rho[i]]*(T+1),
								'k_r':			[k_r[i]]*(T+1),
								'lamda':		[lamda[i]]*(T+1),
								'k_l': 			[k_l[i]]*(T+1),
								'eta': 			[eta[i]]*(T+1),
								'alpha_B': 		[alpha_B[i]]*(T+1),
								'B_max': 		[B_max[i]]*(T+1),
								'k_alpha': 		[k_alpha[i]]*(T+1),
								'delta_B': 		[delta_B[i]]*(T+1),
								'delta_Bext': 	[delta_Bext[i]]*(T+1),
								'epsilon_B': 	[epsilon_B[i]]*(T+1),
								'xi_B': 		[xi_B[i]]*(T+1),
								'delta_A': 		[delta_A]*(T+1),
								'epsilon_A':	[epsilon_A[i]]*(T+1),
								'xi_A': 		[xi_A[i]]*(T+1),
								'Vmax_B': 		[Vmax_B[i]]*(T+1),
								'k_B': 			[k_B[i]]*(T+1),
								'Vmax_Bext': 	[Vmax_Bext[i]]*(T+1),
								'k_Bext': 		[k_Bext[i]]*(T+1),
								'Vol': 			[Vol[i]]*(T+1),
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

		fig, ax = pyplot.subplots(2,2, figsize=(16,9))

		t_series = self.data[self.data['strain']==self.data['strain'].unique()[0]]['t'].values
		ax[0,0].plot(t_series, [1e6]*len(t_series), linestyle='-', color='lightgray')
		ax[1,0].plot(t_series, [0.0]*len(t_series), linestyle='-', color='lightgray')

		strains 		= self.data['strain'].unique()
		strain_palette 	= seaborn.color_palette("husl", len(strains))
		for i, strain in enumerate(strains):

			t_series 	= self.data[self.data['strain']==strain]['t'].values
			N_series 	= self.data[self.data['strain']==strain]['N'].values
			r_series 	= self.data[self.data['strain']==strain]['r'].values
			l_series 	= self.data[self.data['strain']==strain]['l'].values
			Aext_series = self.data[self.data['strain']==strain]['Aext'].values
			A_series 	= self.data[self.data['strain']==strain]['A'].values
			Bext_series = self.data[self.data['strain']==strain]['Bext'].values
			B_series 	= self.data[self.data['strain']==strain]['B'].values
			k_l_series 	= self.data[self.data['strain']==strain]['k_l'].values
			mediaVolume_series 	= self.data[self.data['strain']==strain]['media_volume'].values
			
			
			if(plot_density):
				ax[0,0].plot(t_series, N_series/mediaVolume_series, linestyle='-', color=strain_palette[i])
				ax[0,0].set_ylabel("population density (cells/mL)")
			else:
				ax[0,0].plot(t_series, N_series, linestyle='-', color=strain_palette[i])
				ax[0,0].set_ylabel("population size (num cells)")
			
			ax[1,0].plot(t_series, r_series, linestyle='--', color=strain_palette[i], alpha=0.5)
			ax[1,0].plot(t_series, l_series, linestyle=':', color=strain_palette[i], alpha=0.5)
			ax[1,0].plot(t_series, r_series-l_series, linestyle='-', color=strain_palette[i])
			ax[1,0].set_ylabel("growth and lysis rates")

			ax[0,1].plot(t_series, k_l_series, linestyle=':', color=strain_palette[i], alpha=0.5)
			ax[0,1].plot(t_series, A_series, linestyle='-', color=strain_palette[i])
			ax[0,1].set_ylabel("antibiotic concentration (ug/mL)")
			
			ax[1,1].plot(t_series, B_series, linestyle='-', color=strain_palette[i])
			ax[1,1].plot(t_series, Bext_series, linestyle='--', color=strain_palette[i])
			ax[1,1].set_ylabel("B-lactamase concentration")
			
		ax[0,1].plot(t_series, Aext_series, linestyle='--', color='black')

		seaborn.set_style("ticks")
		seaborn.despine()

		fig.tight_layout()

		pyplot.show()

