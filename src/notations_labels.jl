
# generating symbols 
symb_si(s,i) = Symbol(s*string(i))
symb_Δm(i) = symb_si("Δm",i)
symb_m(i) = symb_si("m",i)
symb_Γ(i) = symb_si("Γ",i)

L_ΞcKπ = L"m(\varXi_c^{+}K^-\pi^-)\,\,\,[\mathrm{GeV}]"
L_ΞcK = L"m(\varXi_c^{+}K^-)\,\,\,[\mathrm{GeV}]"
L_Kπ = L"m(K^-\pi^-)\,\,\,[\mathrm{GeV}]"
L_ΞcK_Ξc_K = L"m(\varXi_c^{+}K^-) - m_{\varXi_c^{+}} - m_{K^-}\,\,\,[\mathrm{MeV}]"
L_cand = L"\mathrm{Number\,\,of\,\,candidates}"
L_cosθ = L"\cos\,\theta"

L_cand_ofXeV(bin_width, units="MeV") = latexstring("\\mathrm{Candidates}\\,/\\,$(bin_width)\\,\\mathrm{$(units)}")
L_Ωctitle(n) = latexstring("\\varOmega_c($(n))^0")
