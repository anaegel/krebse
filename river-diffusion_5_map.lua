--------------------------------------------------------------------------------
-- Example of a problem on 1D line segments in 2D
-- (C) G-CSC, Uni Frankfurt, 2021
--------------------------------------------------------------------------------
ug_load_script("ug_util.lua")
ug_load_script("util/profiler_util.lua")
ug_load_script("util/load_balancing_util.lua")

local gridName = util.GetParam("-grid", "grids/simple-river-y2-meters.ugx")
-- local gridName = util.GetParam("-grid", "grids/map_scale_1km_cut_2_wire_cleanup_meters.ugx")

local charLength = 3.0

local ARGS = {}
ARGS.eps = 1.0 -- diffusion constant
ARGS.dim = util.GetParamNumber("-dim", 2, "dimension")
ARGS.numPreRefs = util.GetParamNumber("-numPreRefs", 0, "number of refinements before parallel distribution")
ARGS.numRefs    = util.GetParamNumber("-numRefs",    5, "number of refinements") 
ARGS.startTime  = util.GetParamNumber("-start", 0.0, "start time")
ARGS.endTime    = util.GetParamNumber("-end", 0.2*charLength*charLength/ARGS.eps, "end time")
ARGS.endTime    = util.GetParamNumber("-end", 365, "end time")
ARGS.dt 	= util.GetParamNumber("-dt", ARGS.endTime*5e-4, "time step size")
util.CheckAndPrintHelp("Time-dependent problem setup example");
print(" Choosen Parameter:")
print("    numRefs      = " .. ARGS.numRefs)
print("    numPreRefs   = " .. ARGS.numPreRefs)
print("    startTime 	= " .. ARGS.startTime)
print("    endTime 	= " .. ARGS.endTime)
print("    dt 		= " .. ARGS.dt)
print("    grid         = " .. gridName)

-- choose algebra
InitUG(ARGS.dim, AlgebraType("CPU", 3));

-- Create, Load, Refine and Distribute Domain
local mandatorySubsets = {"River", "Sink", "Source"}
local dom = nil
if withRedist == true then
	dom = util.CreateDomain(gridName, 0, mandatorySubsets)
	balancer.RefineAndRebalanceDomain(dom, ARGS.numRefs)
else
	dom = util.CreateAndDistributeDomain(gridName, ARGS.numRefs, ARGS.numPreRefs, mandatorySubsets)
end
print("\ndomain-info:")
print(dom:domain_info():to_string())

-- Labels
require "$labels"
local label = LABELS

-- create Approximation Space
print(">> Create ApproximationSpace")
local approxSpace = ApproximationSpace(dom)
approxSpace:add_fct(label.inv_A, "Lagrange", 1)
approxSpace:add_fct(label.sus_B, "Lagrange", 1)
approxSpace:add_fct(label.inf_B, "Lagrange", 1)

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
-- local velocity = ScaleAddLinkerVector()
-- velocity:add(0.1,downStreamVec)
local Amax = 2.0
local Bmax = 2.0

-- Species A - invasive
-- Species B - native

-- My variables
require "$params"
local vars = PARAMS

local upwind = FullUpwind()
local elemDisc = {}
elemDisc[label.inv_A] = ConvectionDiffusionFV1(label.inv_A, "River")
elemDisc[label.sus_B] = ConvectionDiffusionFV1(label.sus_B, "River")
elemDisc[label.inf_B] = ConvectionDiffusionFV1(label.inf_B, "River")

-- Generelle Variablen für die Gleichung von Species A
-- local growth_a = vars.r_a * (1.0 - (vars.k_a * elemDisc[label.inv_A]:value() + vars.k_b * (elemDisc[label.sus_B]:value() + elemDisc[label.inf_B]:value()))/vars.k_max);
-- local death_a = vars.m_a;
local growth_a = vars.r_a * (1.0 - (vars.k_a * elemDisc[label.inv_A]:value()))

-- Gleichung für invasive A
elemDisc[label.inv_A]:set_mass_scale(1.0)      -- \partial A / \partial t
elemDisc[label.inv_A]:set_diffusion(vars.d_a)  -- \partial_x (d_a* \partial_x A)
elemDisc[label.inv_A]:set_reaction_rate(growth_a)
-- elemDisc[label.inv_A]:set_reaction_rate(elemDisc["A"]:value()-Amax) -- Vermehrung: (Amax-A)*A

-- Generelle Variablen für die Gleichung von Species B
local growth_b = vars.r_b * (1.0 - (vars.k_a * elemDisc[label.inv_A]:value() + vars.k_b * (elemDisc[label.sus_B]:value() + elemDisc[label.inf_B]:value()))/vars.k_max);
local death_b = vars.m_b;
local infect_b = vars.i_ab * elemDisc[label.inv_A]:value() + vars.i_bb * elemDisc[label.inf_B]:value();

-- Gleichung fuer susceptable B
elemDisc[label.sus_B]:set_mass_scale(1.0) 
elemDisc[label.sus_B]:set_diffusion(vars.d_b)
-- elemDisc[label.sus_B]:set_upwind(upwind)
elemDisc[label.sus_B]:set_reaction_rate(growth_b - infect_b)  	-- Sterben: -B*A 

-- Gleichung für infected B
elemDisc[label.inf_B]:set_mass_scale(1.0) 
elemDisc[label.inf_B]:set_diffusion(vars.d_b)
elemDisc[label.inf_B]:set_reaction_rate(0.0 - death_b)
-- elemDisc[label.inf_B]:set_reaction((0.0 - infect_b) * elemDisc[label.sus_B]:value())
-- elemDisc[label.inf_B]:set_reaction_rate(infect_b/elemDisc[label.inf_B]:value() - death_b)  	-- Sterben: -B*A 

-- Boundary Conditions
require "$boundary"
local dirichletBND = BOUNDARY
--dirichletBND:add(Amax, label.inv_A, "Sink")	-- An der Muendung nur invasive Spezies A
--dirichletBND:add(Bmax, label.sus_B, "Source")	-- An Quellen wird Spezies B ausgesetzt.
--dirichletBND:add(Bmax, label.inf_B, "Sink") -- An Muendung wird Spezies B ausgesetzt.

-- Create discretization.
local domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc[label.inv_A])
domainDisc:add(elemDisc[label.sus_B])
domainDisc:add(elemDisc[label.inf_B])
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
    linSolver = LU(), --[[{
        type = "bicgstab",
        precond = {
        type = "ilu",
        -- sort = true,
        },
        convCheck = {
            type	= "standard",
            iterations	= 100,
            absolute	= 1e-11,
            reduction	= 1e-13,
            verbose=true
        },
	},--]]
    lineSearch = {
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
Interpolate(0.0, u, label.inv_A, ARGS.startTime)   		    -- no predators
Interpolate(Bmax*0.0, u, label.sus_B, ARGS.startTime)  		-- only prey
-- TODO: Do I need one for inf_B?

-- Configure VTK output.
local vtk = VTKOutput()
local vtkObserver=VTKOutputObserver("zMyFile.vtk", vtk)
-- local vtkObserver=VTKOutputObserver("MyFile.vtk", vtk, 3) -- Zeiiteinheiten

vtk:select(label.inv_A, "Species A - Invasive")
vtk:select(label.sus_B, "Species B - Susceptable")
vtk:select(label.inf_B, "Species B - Infected")

local my_count = 2
function outputvtk(u, step, time)
  if step % my_count == 0 then
    vtk:print("vtk/krebse", u, step, time)
  end
end

-- Perform time stepping loop.
-- util.SolveNonlinearTimeProblem(u, domainDisc, solver, outputvtk, "vtk/Krebse",
--						   "ImplEuler", 1, ARGS.startTime, ARGS.endTime, ARGS.dt, ARGS.dt * 1e-6);


local limexLSolver = nil
local limexNLSolver = nil

local limexConvCheck=ConvCheck(1, 5e-8, 1e-10, true)
limexConvCheck:set_supress_unsuccessful(true) 
 
--[[ if false then 
 limexLSolver = {}
 limexNLSolver = {}
for i=1,nstages do 

  limexLSolver[i] = util.solver.CreateSolver(solverDesc)
    
  limexNLSolver[i] = NewtonSolver()
  limexNLSolver[i]:set_linear_solver(limexLSolver[i])
  limexNLSolver[i]:set_convergence_check(limexConvCheck)
  
  print(limexNLSolver[i])
end
else --]]
limexLSolver = util.solver.CreateSolver(solverDesc.linSolver)
limexNLSolver = NewtonSolver()
 limexNLSolver:set_linear_solver(limexLSolver)
 limexNLSolver:set_convergence_check(limexConvCheck)
  print(limexNLSolver)
--[[ end --]]


-- local refObserver = PlotRefOutputObserver("DirichletValue", vtk) -- now obsolete
local luaObserver = LuaCallbackObserver()

-- work-around (waiting for implementation of SmartPtr forward to lua...)
function luaPostProcess(step, time, currdt)
  print("LUAPostProcess: "..step..","..time..","..currdt)
  -- TODO: vtk output
  -- postProcess(luaObserver:get_current_solution(), step, time, currdt)
  return 0;
end
luaObserver:set_callback("luaPostProcess")

--  Euclidean (algebraic) norm
--local estimator = Norm2Estimator()  
--tol = 0.37/(gridSize)*tol


--print (estimator)
local limexEstimator = CompositeGridFunctionEstimator()
--limexEstimator:add(L2ComponentSpace(label.inf_B, 2))  -- Genauigkeit der Quadraturformel,, norm im finite element raum
--limexEstimator:add(L2ComponentSpace(label.sus_B, 2))  
limexEstimator:add(H1SemiComponentSpace(label.inv_A, 2))  

-- descriptor for integrator
local limexDesc = {

  nstages = 2, -- stufen in neville schema
  steps = {1,2,3,4,5,6}, -- unterteilungen
  domainDisc=domainDisc,
  nonlinSolver = limexNLSolver,
  
  tol = 1e-3, -- toleranz 1e-2
  dt = ARGS.dt, --dtlimex,
  dtmin = ARGS.dt * 1e-6, --- minimale zeitschrittweite
  
  rhoSafetyOPT = 0.8, -- sicherheitsfaktor, kleiner Eins
  
}


-- setup for time integrator
local limex = util.limex.CreateIntegrator(limexDesc)

--limex:set_dt_min(1e-9)
limex:add_error_estimator(limexEstimator)
limex:set_increase_factor(2.0) -- maximale erhöhung, von schrittweite in einem step

if (vtk) then 
   limex:attach_observer(vtkObserver)
end

limex:attach_observer(luaObserver)
-- limex:attach_observer(refObserver)

-- limex:set_stepsize_greedy_order_factor(1.0)
-- limex:select_cost_strategy(LimexNonlinearCost())
-- limex:disable_matrix_cache()  -- recompute ()
limex:enable_matrix_cache() -- keep matrix 

-- solve problem

-- print(">> Peclet number:"..50.0*1.0/eps)
-- print(">> Grid Peclet number:"..50.0*gridSize/eps)
-- print(">> Solve problem")
local cstart=os.clock()
limex:apply(u, ARGS.endTime, u, ARGS.startTime)
local cend=os.clock()
-- print ("CDELTA=\t"..cend - cstart)

 
print("Writing profile data")
WriteProfileData("profile_data.pdxml")
util.PrintProfile_TotalTime("main ")
-- end group app_convdiff
--[[!
\}
]]--
