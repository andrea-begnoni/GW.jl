using Documenter
using Pkg
using ForwardDiff

current_dir = @__DIR__

# Get the parent directory
parent_dir = dirname(current_dir)

# Print the parent directory

Pkg.activate(parent_dir) # activate the environment in the current directory
# this is needed to use this package without the need to install it from the registry (where the package is not yet registered)
using GW

makedocs(sitename="My Documentation")