abstract type SuperType

end

function SuperTypeMethod1(this::SuperType)
    return this.super_field
end

mutable struct Subtype1 <: SuperType
    super_field
end

mutable struct Subtype2 <: SuperType
    
end

# The abstract type is purely for dynamic dispatch, to share the field from the 
# supertype, all subtype must share the same internal data. 