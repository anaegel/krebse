require "$labels"

local label = LABELS

-- My variables
BOUNDARY = DirichletBoundary()
BOUNDARY:add(0.5, label.inv_A, "Sink") -- An Muendung wird Spezies A ausgesetzt.

