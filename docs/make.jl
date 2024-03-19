push!(LOAD_PATH,"../src/")
using ASCertain
using Documenter

makedocs(
         sitename = "ASCertain.jl",
         modules  = [ASCertain],
         pages=[
                "Home" => "index.md"
               ])

deploydocs(;
           repo="github.com/darnstrom/ASCertain.jl",
          )
