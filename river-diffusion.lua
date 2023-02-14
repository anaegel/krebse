--------------------------------------------------------------------------------
-- Example of a problem on 1D line segments in 2D
-- (C) G-CSC, Uni Frankfurt, 2021
--------------------------------------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("util/profiler_util.lua")
ug_load_script("util/load_balancing_util.lua")


--local gridName = util.GetParam("-grid", "grids/simple-river-y.ugx")
local gridName = util.GetParam("-grid", "grids/simple-river-y.ugx")
local charLength = 3.0


local ARGS = {}
ARGS.eps = 1.0 -- diffusion constant

ARGS.dim = util.GetParamNumber("-dim", 2, "dimension")
ARGS.numPreRefs = util.GetParamNumber("-numPreRefs", 0, "number of refinements before parallel distribution")
ARGS.numRefs    = util.GetParamNumber("-numRefs",    6, "number of refinements") 
ARGS.startTime  = util.GetParamNumber("-start", 0.0, "start time")
ARGS.endTime    = util.GetParamNumber("-end", 0.2*charLength*charLength/ARGS.eps, "end time")
ARGS.dt 		   = util.GetParamNumber("-dt", ARGS.endTime*5e-4, "time step size")

util.CheckAndPrintHelp("Time-dependent problem setup example");


print(" Choosen Parameter:")
print("    numRefs      = " .. ARGS.numRefs)
print("    numPreRefs   = " .. ARGS.numPreRefs)
print("    startTime 	= " .. ARGS.startTime)
print("    endTime 		= " .. ARGS.endTime)
print("    dt 			= " .. ARGS.dt)
print("    grid         = " .. gridName)

-- choose algebra
InitUG(ARGS.dim, AlgebraType("CPU", 2));

-- Create, Load, Refine and Distribute Domain
local mandatorySubsets = {"River", "Sink"}
local dom = nil
if withRedist == true then
	dom = util.CreateDomain(gridName, 0, mandatorySubsets)
	balancer.RefineAndRebalanceDomain(dom, ARGS.numRefs)
else
	dom = util.CreateAndDistributeDomain(gridName, ARGS.numRefs, ARGS.numPreRefs, mandatorySubsets)
end

print("\ndomain-info:")
print(dom:domain_info():to_string())

-- create Approximation Space
print(">> Create ApproximationSpace")
local approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("A", "Lagrange", 1)
approxSpace:add_fct("B", "Lagrange", 1)

-- Order downstream
local downStreamVec = EdgeOrientation(dom)
OrderDownwind(approxSpace, downStreamVec, 3.14/20.0)

-- OrderLex(approxSpace, "x");			-- Lexicographic order of indices. 

--------------------------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
--------------------------------------------------------------------------------

print (">> Setting up Assembling")

--------------------------------------------------------------------------------
-- Gleichungen:
--------------------------------------------------------------------------------
--  \partial_t a + \partial_x [-lambda  \partial_x a + v a] + (amax-a)*a = 0 
--  \partial_t b + \partial_x [-lambda  \partial_x b + v b] + a*b = 0

local lambda  = 1.0
local velocity = ScaleAddLinkerVector()
velocity:add(0.1,downStreamVec)

local Amax = 2.0
local Bmax = 1.0

local upwind = FullUpwind()
local elemDisc = {}
elemDisc["A"] = ConvectionDiffusionFV1("A", "River")
elemDisc["B"] = ConvectionDiffusionFV1("B", "River")


-- Gleichung für A
elemDisc["A"]:set_mass_scale(1.0)                   		-- \partial A / \partial t
elemDisc["A"]:set_diffusion(1.0)                        -- \partial_x (lambda* \partial_x A)
elemDisc["A"]:set_velocity( downStreamVec)   						-- \partial_x (v*A)
elemDisc["A"]:set_upwind(upwind)
elemDisc["A"]:set_reaction_rate(elemDisc["A"]:value()-Amax) -- Vermehrung: (Amax-A)*A

-- Gleichung fuer B
elemDisc["B"]:set_mass_scale(1.0) 
elemDisc["B"]:set_diffusion(lambda)                      	-- \partial_x (lambda* \partial_x B)
elemDisc["B"]:set_velocity(downStreamVec)      						-- \partial_x (v*B)
elemDisc["B"]:set_upwind(upwind)
elemDisc["B"]:set_reaction_rate(elemDisc["A"]:value())  	-- Sterben: -B*A 


-- 
local dirichletBND = DirichletBoundary()
dirichletBND:add(Amax, "A", "Sink")						-- An der Muendung nur invasive Spezies A
dirichletBND:add(Bmax/10, "B", "Source2")			-- An Quelle 2 wird Spezies B ausgesetzt.


-- Create discretization.
local domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc["A"])
domainDisc:add(elemDisc["B"])
domainDisc:add(dirichletBND)



--------------------------------------------------------------------------------
--  Algebra
--------------------------------------------------------------------------------
print (">> Setting up Algebra Solver")

-- if the non-linear problem shall should be solved, one has to wrap the solver
-- inside a newton-solver. See 'solver_util.lua' for more details.
solverDesc = 
{
  type = "newton",
  
 
  linSolver = 
  {
	 type = "bicgstab",
	 
	  precond = {
	   type = "ilu",
	   -- sort = true,
		},
	
		convCheck = {
			type		= "standard",
			iterations	= 100,
			absolute	= 1e-11,
			reduction	= 1e-13,
			verbose=true
		},
	
	},
	
  lineSearch =
  {
    type = "standard",
      maxSteps    = 10,
      lambdaStart   = 1.0,
      lambdaReduce  = 0.5,
      acceptBest    = false,
      checkAll    = false
  },	

}

local dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_vtk_output(false)

local solver = util.solver.CreateSolver(solverDesc)
--solver:set_debug(dbgWriter)

--------------------------------------------------------------------------------
--  Apply Solver
--------------------------------------------------------------------------------

-- Set initial value.
print(">> Interpolation start values")
local u = GridFunction(approxSpace)
Interpolate(0.0, u, "A", ARGS.startTime)   				-- no predators
Interpolate(Bmax*0.0, u, "B", ARGS.startTime)  		-- only prey


-- Configure VTK output.
local vtk = VTKOutput() 
vtk:select("A", "Species A")
vtk:select("B", "Species B")

-- Perform time stepping loop.
util.SolveNonlinearTimeProblem(u, domainDisc, solver, vtk , "vtk/Krebse",
							   "ImplEuler", 1, ARGS.startTime, ARGS.endTime, ARGS.dt, ARGS.dt * 1e-6); 


print("Writing profile data")
WriteProfileData("profile_data.pdxml")
util.PrintProfile_TotalTime("main ")

-- end group app_convdiff
--[[!
\}
]]--
