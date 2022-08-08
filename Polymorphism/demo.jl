### We demonstrate emulating OOP in julia via mixers package and dynamic dispatch types hierachy. 

using Mixers

abstract type AbstractAnimal

end

@premix mutable struct AbstractAnimalTemplate
    domain::String
    kingdom::String
    phylum::String
    class::String
    order::String
    family::String
    genus::String
    species::String

end

"""
    This method should be inherited for other subtypes by using the correct call signature for the second parameter. 
"""
function GenerateSummary(this::AbstractAnimal, ::Type{AbstractAnimal})
    summary = Vector{String}()
    push!(summary, "Taxonomic Ranks: ====== \n")
    for sym in [:domain, :kingdom, :phylum, :class, :order, :family, :genus, :species]
        push!(
            summary, "$(sym) = $(getfield(this, sym))\n"
        )
    end
    push!(summary, "END Taxonomic Ranks ======")
return foldl(*, summary) end


@AbstractAnimalTemplate mutable struct HomoSapient <: AbstractAnimal
    name::String
    jobs::Vector{String}
    education_highest::Union{Nothing, String}
    

    function HomoSapient()
        this = new()
        this.domain = "eukaryota"
        this.kingdom = "animalia"
        this.phylum = "chordata"
        this.class = "mammalia"
        this.order = "primates"
        this.family = "hominidea"
        this.genus = "homo"
        this.species = "sapient"

        this.name = "Homer Simpson"
        this.jobs = ["Spring Field's Comedian", "Bart's father", "Marge's Husband"]
        this.education_highest = "Springfield University, Degree in Nuclear Physics"

    return this end
end

"""
    This method can now inherit from abstract methods of its own, or other supertypes by the correct type 
    singature. 
"""
function GenerateSummary(this::HomoSapient)
    result = GenerateSummary(this, AbstractAnimal)
    result *= "\nThey think they are something important but they are really not.\n "
    result *= "Individual Profile ======\n"
    result *= "name:"*this.name*"\n"
    result *= "jobs:"*foldl(*, this.jobs)*"\n"
    result *= "education_highest: "*this.education_highest
return result end


@AbstractAnimalTemplate mutable struct RedFox
    name::String
    location::String
    function RedFox()
        this = new()
        # Try simplifying this part with polymorphism too. 
        this.domain = "eukaryota"
        this.kingdom = "animalia"
        this.phylum = "chordata"
        this.class = "mammalia"
        this.order = "canivora"
        this.family = "canidea"
        this.genus = "vulpes"
        this.species = "vulpes vulpes"

    return this end
end

function GenerateSummary(this::RedFox)
    result = GenerateSummary(this, AbstractAnimal)
    result *= "\n They are cute animals in homo sapeint's views."
return result end



homer_simpson = HomoSapient()
GenerateSummary(homer_simpson)|>print