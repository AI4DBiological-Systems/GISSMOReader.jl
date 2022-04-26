module GISSMOReader

import Unicode
import JLD
import InfoZIP
import DelimitedFiles

import HTTP
import JSON3
import JSON

include("../src/database/download_edata.jl")
include("../src/database/helpers.jl")
include("../src/database/BMRB_rest.jl")
include("../src/database/GISSMO_entries.jl")
include("../src/database/metadata.jl")

include("../src/assemble.jl")
include("../src/utils.jl")

end
