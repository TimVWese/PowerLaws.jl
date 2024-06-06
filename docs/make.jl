using Documenter
push!(LOAD_PATH,"../src/")
using PowerLaws

makedocs(
    sitename = "PowerLaws",
    format = Documenter.HTML(),
    modules = [PowerLaws],
    pages = [
        "Home" => "index.md",
        "Reference" => "reference.md"
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
