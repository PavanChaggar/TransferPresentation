### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# ╔═╡ 4cf73862-25c7-11ec-174b-6dcf37428913
begin
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ 5a0591b7-1a34-4e3e-9862-c772fc3159f4
html"""<style>
main {
max-width: 900px;
}"""

# ╔═╡ 18225564-8512-4fca-87c8-a95ec2fa0d05
html"<button onclick='present()'>present</button>"

# ╔═╡ 4ff67e50-ccdd-479f-8280-e04ab2354ce4
md" 
# Transfer Thesis Viva 

**Pavanjit Chaggar, October 2021**

Supervised by Alain Goriely and Saad Jbabdi, with support from Stefano Magon and Gregory Klein at Roche

Presentation: https://github.com/PavanChaggar/TransferPresentation"

# ╔═╡ abc58f7f-c4c1-47b6-861a-ab679d34bc95
md" 
# Overview and Introduction

- Alzheimer's disease (AD)
- Mathematical models of AD
- Inference Workflow
- Preliminary results
"

# ╔═╡ 95d6223a-c12e-4b26-8c4e-d59a59c7d129
md" 
# Alzheimer's Disease -- A Brief Summary
Alzheimer's is characterised by gradual neurodegeneration associated with pervasive spreading of toxic protein species. 
In particular, two proteins, Amyloid beta (Aβ) and tau-protein (τP) are believed to underlie and drive the development of pathology. 
Historically, Aβ was primarily invstigated as the primary cause of AD. However More recent work has focussed on τP, in part because it spreads very predictably and is more tightly coupled with atrophy and symotom onset.
"

# ╔═╡ d75eb4e7-2fbf-44ca-af86-bf67fc1d393d
md" 
## A Pernicious Pair of Predictable Prion Proteins
Both Aβ and τP grow via an autocatalytic process resembling those displayed by prions. 
This process is summarised as: 
"

# ╔═╡ a0cb7614-2ab3-44d1-9202-02f19915edf6
html"""
<img src="https://github.com/PavanChaggar/inference-methods-UCSF0721/blob/main/assets/images/heterodimerkinetics.png?raw=true" height=250 width=500 vspace=50, hspace=175>"""

# ╔═╡ f45c2cd6-baf6-4ce3-84b0-5bf8fb9e67d4
md"## Braak Stages of Tau protein
In most AD cases, τP follows a predictable pattern of spreading, starting in the entorhinal cortex before spreading through the hippocampal regions, lateral cortex and finally into the neocortex. Atrophy tends to closely follow the spreading pattern of Tau, more so than that of Aβ"


# ╔═╡ Cell order:
# ╠═4cf73862-25c7-11ec-174b-6dcf37428913
# ╠═5a0591b7-1a34-4e3e-9862-c772fc3159f4
# ╠═18225564-8512-4fca-87c8-a95ec2fa0d05
# ╟─4ff67e50-ccdd-479f-8280-e04ab2354ce4
# ╟─abc58f7f-c4c1-47b6-861a-ab679d34bc95
# ╟─95d6223a-c12e-4b26-8c4e-d59a59c7d129
# ╟─d75eb4e7-2fbf-44ca-af86-bf67fc1d393d
# ╟─a0cb7614-2ab3-44d1-9202-02f19915edf6
# ╠═f45c2cd6-baf6-4ce3-84b0-5bf8fb9e67d4
