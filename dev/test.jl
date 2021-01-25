(var"##MTIIPVar#256", var"##MTKArg#254")->begin
        @inbounds begin
                begin
                    (ModelingToolkit.fill_array_with_zero!)(var"##MTIIPVar#256")
                    let (x, y) = (var"##MTKArg#254"[1], var"##MTKArg#254"[2])
                        var"##MTIIPVar#256"[1] = (+)((^)(x, 2), y)
                        var"##MTIIPVar#256"[2] = (+)((^)(y, 2), x)
                    end
                end
            end
        nothing
    end