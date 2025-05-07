
function pretty_print(hp::HyperParameters)
    println("HyperParameters:")
    println("================")
    for field in fieldnames(HyperParameters)
        val = getfield(hp, field)
        # Format vectors nicely
        val_str = isa(val, AbstractVector) ? "[" * join(val, ", ") * "]" : string(val)
        @printf("  %-22s : %s\n", string(field), val_str)
    end
end
