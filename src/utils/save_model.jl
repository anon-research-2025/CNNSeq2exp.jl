function model2cpu(m::model)
    img_filters_cpu = 
        [learned_codeimg_filters(
            m.img_filters[i].f |> Array, m.img_filters[i].eta)
            for i in axes(m.img_filters, 1)]
    pwms_cpu = learned_pwms(
                    m.pwms.f |> Array, 
                    m.pwms.b |> Array, 
                    m.pwms.w |> Array,
                    m.pwms.etas)

    return model(pwms_cpu, 
                 img_filters_cpu, 
                 m.w |> Array, 
                 m.beta)
end

function model2gpu(m::model)
    img_filters_gpu = 
        [learned_codeimg_filters(
            m.img_filters[i].f |> gpu, m.img_filters[i].eta)
            for i in axes(m.img_filters, 1)]
    pwms_gpu = learned_pwms(
                    m.pwms.f |> gpu, 
                    m.pwms.b |> gpu, 
                    m.pwms.w |> gpu,
                    m.pwms.etas)

    return model(pwms_gpu, 
                 img_filters_gpu, 
                 m.w |> gpu, 
                 m.beta)
end