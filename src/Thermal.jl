## THERMAL LOAD MODELLING
# SIMPLIFIED METHOD BY UNI EN ISO 13790:2008  UNI EN ISO 52016-1:2017

using DifferentialEquations, LinearAlgebra

## STANDARD AND FIXED VALUES

    # CAPACITY per unit area by building element by regulation

"""
Inserire modalità di dizionario che per ogni tipo di elemento prende il giusto valore di K indicato dalla norma, e vedere i valori della norma
FARE UN ELENCO O DIZIONARIO, NON SO, IN CUI SCELGO IL NUMERO DEGLI ELEMENTI OPACHI “OPAQUE” O IL NUMERO DEGLI ELEMENTI FINESTRATI “WINDOW” DALL’ELENDO DEL NUEMRO DEGLI ELEMENTI PRESENTI IN OGNI EDIFICIO. AD OGNI ELEMENTO SPECIFICO ASSOCIARE LA SUA AREA, LA SUA CAPACITA’ INTERNA PER UNITA’ DI SUPERFICIE, IL SUO COEFFICIENTE DI SCAMBIO TERMICO
"""

# In ECModels.jl
"""
load data 
ECModel.data = raw_data["data"]
ECModel.gen_data = general(ECModel.data)
ECModel.market_data = market(ECModel.data)
ECModel.users_data = users(ECModel.data)
"""

# In utils.jl

"""
    @enum BUILDING_CLASS

Enumeration type to specify the class of the buildings.
Implemented values:
- VLIGHT: very light class
- LIGHT: light class
- MED: medium class
- HEAVY: heavy class
- VHEAVY: very heavy class
"""
@enum BUILDING_CLASS VLIGHT=0 LIGHT=1 MED=2 HEAVY=3 VHEAVY=4

ALL = collect(instances(BUILDING_CLASS))  # all classes code

# DEVICES = setdiff(ALL, [LOAD])  # devices codes
# GENS = [REN, THER]  # generator codes

class_codes = Base.Dict("very_light"=>VLIGHT, "light"=>LIGHT,"medium"=>MED,"heavy"=>HEAVY, "very_heavy"=>VHEAVY)




"Function to safely get a field of a dictionary with default value"
@inline field_d(d::AbstractDict, field, default=nothing) = (field in keys(d) ? d[field] : default)
@inline field_i(d, field) = field_d(d, field, 0)
@inline field_f(d, field) = field_d(d, field, 0.0)
"Function get field that throws an error if the field is not found"
@inline function field(d::AbstractDict, field, desc=nothing)
    if field in keys(d)
        return d[field]
    else
        msg = isnothing(desc) ? "Field $field not found in dictionary $(keys(d))" : desc
        throw(KeyError(msg))
    end
end

"Function to get the general parameters"
general(d::AbstractDict) = field(d, "general")
"Function to get the users configuration"
users(d::AbstractDict) = field(d, "users")
"Function to get the market configuration"
market(d::AbstractDict) = field(d, "market")

"Function to get the profile dictionary"
profiles(d::AbstractDict) = field_d(d, "profile")

"Auxiliary function to check if the key 'class' is available in the dictionary d, otherwise false"
has_class(d::AbstractDict) = ("class" in keys(d))
has_class(d) = false  # if d is not an abstract dictionary, then return false

"Function to get the components list of a dictionary"
function property(d::AbstractDict)
    return Dict(k=>v for (k,v) in d if has_class(v))
end
"Function to get the property value of a dictionary"
property(d, p_name) = field(property(d), p_name)
"Function to get the property value of a dictionary"
field_property(d, p_name, pf_name) = field(property(d, p_name), pf_name)
"Function to get the properties value of a dictionary, with default value"
field_property(d, p_name, pf_name, default) = field_d(component(d, c_name), pf_name, default)
"Function to know if a dictionary has a particular property"
has_property(d, p_name, pf_name) = haskey(property(d, p_name), pf_name)

"Function to get a specific profile"
function profile(d, profile_name)
    profile_block = profiles(d)
    return field(profile_block, profile_name)
end
"Function to get a specific profile"
function profile_component(d, c_name, profile_name)
    profile_block = profiles(component(d, c_name))
    return field(profile_block, profile_name)
end



"Function to get the building class of a property"
building_class(d, prop_name) = class_codes[field(property(d, prop_name), "class")]

"Function to get the list of the buildings, referred to a specific class of building, for a user"
function building_names(d, b_class::BUILDING_CLASS)
    props = property(d)
    return [bc for bc in keys(props) 
        if asset_type(comps, bc) == b_class]
end

"Function to get the list of the all the buildings for a user"
building_names(d) = collect(keys(property(d)))

"Function to get the list of the properties, referred to a -one or more- particoular building, for a user, as a list of elements"
function building_names(d, b_class::Vector{BUILDING_CLASS})
    props = property(d)
    return [bc for bc in keys(props) 
        if building_class(comps, bc) ∈ b_classes]
end


"A distinction could be made between public and private buildings"
" replace DEVICES with PUBLIC and GENS with PRIVATE "

"Function to get the list of devices for a user" -> "Function to get the list of public buildings for a user"
device_names(d) = asset_names(d, DEVICES)

"Function to get the list of generators for a user" -> "Function to get the list of private buildings for a user"
generator_names(d) = asset_names(d, GENS)

"Function to check whether an user has any asset"
has_any_asset(d, a_types::Vector{ASSET_TYPE}=DEVICES) = !isempty(asset_names(d, a_types))

"Function to check whether an user has an asset type"
has_asset(d, atype::ASSET_TYPE) = !isempty(asset_names(d, atype))

"Function to check whether an user has an asset given its name"
has_asset(d, aname::AbstractString) = aname in keys(d)


"in base_models.jl"

"definitions"
# k # thermal capacity per unit area
# A_m # effective mass area; m2
# U_w # windowed surface trasmission heat transfer coefficient; W/m2K
# A_w # windowed surface area; m2
# h_ms = 9.1 # opaque surface trasmission heat transfer coefficient; W/m2K
# h_is = 3.45 # internal air/internal surface trasmission heat transfer coefficient; W/m2K
# H_op # external air/thermal mass trasmission heat transfer coefficient
# S_at = 4.5 # adimentional ratio between internal surface/heated surface
# A_f # neat heated surface; m2
# rho_air = 1.2 # air density; kg/m3
# c_air = 1005 # air specific heat; J/kgK
# n = 0.4 # air exchange rate; h-1
# Vol_flow_rate_air # air flow rate; m3/h
# Q_int = 3.5 # internal sources, external sources; W/m2
# Q_sol # solar radiation, external sources; W/m2
# T_i # internal temperature, time-depedent; K
# T_int_H_set = 20+273.15 # internal temperature heating set point; K
# T_int_C_set = 26+273.15 # internal temperature cooling set point; K


## Expressions
# DINAMIC ELEMENTS
    # THERMAL CAPACITY by user 
    @expression(model_user, C_m[u in user_set, t in time_set, b in building_class(users_data[u])],
        sum(Float64[field_property(users_data[u], b, "k")*field_property(users_data[u], b, "A")[t] for i in number_element(users_data[u], ELEMENT)])
    )
    # EFFECTIVE MASS AREA by user 
    @expression(model_user, A_m[u in user_set, t in time_set, b in building_class(users_data[u])],
        (C_m^2)/(sum(Float64[field_property(users_data[u], b, "k")^2*field_property(users_data[u], b, "A_m")^2[t] for i in number_element(users_data[u], ELEMENT)]))
        )
#    

# TRASMISSION HEAT TRANSFER COEFFICIENTS
    # WINDOWED SURFACE trasmission heat transfer coefficient by user
    @expression(model_user, H_tr_w[u in user_set, b in building_class(users_data[u])],
        sum(Float64[field_property(users_data[u], b, "U_w")*field_property(users_data[u], b, "A_w") for v in number_window(users_data[u], WINDOW)])
    )
    
    # OPAQUE SURFACE (wall) trasmission heat transfer coefficient 
    
    # THERMAL MASS/INTERNAL SURFACE trasmission heat transfer coefficient by user
    @expression(model_user, H_tr_ms[u in user_set, b in building_class(users_data[u])],
       Float64[h_ms*field_property(users_data[u], b, A_m[u,b])]
    )

    # EXTERNAL AIR/THERMAL MASS trasmission heat transfer coefficient by user
    @expression(model_user, C_m[u in user_set, b in building_class(users_data[u])],
        Float64[1/((1/H_op)-(1/H_tr_ms[u,b]))]
    )
# attenzione definire H_op!!!


# VENTILATION HEAT TRANSFER COEFFICIENTS
 @expression(model_user, H_ve[u in user_set, b in build_type(users_data[u])],
        Float64[[field_property(users_data[u], b, "rho_air")]*[field_property(users_data[u], b, "c_air")]]*sum(Float64[field_property(users_data[u], b, "Vol_flow_rate_air")[t] for v in number_window(users_data[u], WINDOW)])
)

# INTERNAL AIR/INTERNAL SURFACE TRASMISSION HEAT TRANSFER COEFFICIENT
 @expression(model_user, H_tr_is[u in user_set, b in building_class(users_data[u])],
Float64[h_is*A_tot]
)
# ATTENZIONE definire A_tot e capire chi è A_f !!!
# A_tot: total area facing the internal zone
# A_tot = S_at*A_f


    # POWER FLUX: SOLAR RADIATION Q_sol
# DATI ESTERNI, creare lista oraria
    # POWER FLUX: INTERNAL SOURCES Q_int
# DATI ESTERNI, creare lista oraria

    # POWER FLUX through internal air
    @expression(model_user, Q_ia[u in user_set, b in building_class(users_data[u], t in time_set)],
        Float64[Q_int[t]/2]
    )

    # POWER FLUX through thermal mass
    @expression(model_user, Q_m[u in user_set, b in building_class(users_data[u], t in time_set)],
        Float64[(A_m/A_tot)*((Q_int[t]/2)+Q_sol[t])]
    )

    # POWER FLUX through internal surface
    @expression(model_user, Q_st[u in user_set, b in building_class(users_data[u], t in time_set)],
        Float64[(1-(A_m/A_tot)-(H_tr_w[u,b]/(h_ms*A_tot)))*((Q_int[t]/2)+Q_sol[t])]
    )



##  NODAL ANALISY
# LINEAR DINAMIC SYSTEM DEFINITION BY UNI EN ISO 13790:2008  UNI EN ISO 52016-1:2017
# INSERIRE LE 3 EQUAZIONI ??? o meglio inserire direttamente il sistema dinamico con vettori e matrici ???

# Approach 1: Matrix Form
" Matrix System Definition"
M = [(H_tr_ms+H_tr_w)*(H_tr_is+H_ve)+H_tr_is*H_ve]
A = [((((H_tr_is+H_ve)*H_tr_ms*H_tr_ms)/(M))-H_tr_em-H_tr_ms)/C_m]
B = [((H_tr_em+((H_tr_ms*H_tr_w*(H_tr_is+H_ve))/M))/C_m)+(H_tr_ms*H_ve*H_tr_is/M/C_M); 1/C_m; H_tr_ms*(H_tr_is+H_ve)/M/C_m; H_tr_ms*H_tr_is/M/C_m; H_tr_ms*H_tr_is/M/C_m]'
C = [H_tr_ms*(H_tr_is+H_ve)/M; H_tr_ms*H_tr_is/M]
D = [(H_tr_is*H_tr_w/M)+(H_ve*(H_tr_ms+H_tr_is+H_tr_w)/M), (H_tr_is*H_ve/M)+(H_tr_w*(H_ve+H_tr_is)); 0, 0; H_tr_is/M, (H_ve+H_tr_is)/M; (H_tr_ms+H_tr_is+H_tr_w)/M, H_tr_ms/M; (H_tr_ms+H_tr_is+H_tr_w)/M, H_tr_ms/M]'
"Input u(t) Definition"
function u(t)
    return [T_e(t); Q_m(t); Q_st(t); Q_ia; Q_HC(t)]
end
"Differential System Definition"
function state_space!(dx, x, p, t)
    A, B = p  # Extract A,B fro parameters
    dx .= A * x(t) + B * u(t)  # State Equation
end
"Starting Conditions Definition"
#inserire logica in cui, se t è relativo ai mesi invernali, usa l'ingresso x_H_0, altrimenti usa x_C_0, chiamando x_0 l'ingresso scelto
x_H_0 = [20.0]  # Start Condition for Heating-case
x_C_0 = [26.0]  # Start Condition for Cooling-case

# ATTENZIONE!!! DEFINIRE tspan
tspan = (0.0, )  # Time Span
"Creation and Resolution of the Differential System"
prob = ODEProblem(state_space!, x_0, tspan)
sol = solve(prob)
"Output y(t) Definition"
y_vals = [C * sol.u[j] + D * u(sol.t[j]) for j in eachindex(sol.t)]


# Approach 2: Explicit Equations Form
"Dinamic System Definition"
function dynamic_system!(dT_m, T_m, p, t)
    dT_m = (-(H_tr_em+H_tr_ms)*T_m(t) + H_tr_ms*T_s + Q_m(t) - T_e(t)*H_tr_em)/C_m
    (H_tr_ms+H_tr_is+H_tr_w)*T_s - H_tr_ms*T_m(t) - H_tr_is*T_i(t) = Q_st(t) + H_tr_w*T_e(t)
    (H_tr_is+H_ve)*T_i(t) - H_tr_is*T_s - H_ve*T_inj = Q_ia(t) + Q_HC(t)
end
"Starting Conditions Definition"
#inserire logica in cui, se t è relativo ai mesi invernali, usa l'ingresso x_H_0, altrimenti usa x_C_0, chiamando x_0 l'ingresso scelto. usare gradi giorno.
x_H_0 = [20.0]  # Start Condition for Heating-case
x_C_0 = [26.0]  # Start Condition for Cooling-case

# ATTENZIONE!!! DEFINIRE tspan
tspan = (0.0, )  # Time Span
"Creation and Resolution of the Dinamic System"
prob = ODEProblem(dynamic_system!, u0, tspan)
sol = solve(prob)

##  ACTIVATION/DEACTIVATION HEATING SYSTEM LOGIC

# Set HEATING ACTIVATION
for u in user_set
    for b in building_class(users_data[u])
        for t in time_set
# HEATING
            if T_i[u,b,t] <= T_int_H_set[u,b,t]
                set_activation(Q_HC[u,b,t])
# It should be always true that if the internal Temperature is under the Set point value, then a heating flux power should be available
                   # ATTENZIONE definire T_i e T_int_H_set !!!
            end
# SET COOLING ACTIVATION
            if T_i[u,b,t] >= T_int_C_set[u,b,t]
                set_activation(Q_HC[u,b,t])
# It should be always true that if the internal Temperature is upper the Set point value, then a cooling flux power should be available
                   # ATTENZIONE definire T_int_C_set !!!
            end
# SET HEATING/COOLING DEACTIVATION
            if T_i[u,b,t] <= T_int_C_set[u,b,t] && T_i[u,b,t] >= T_int_H_set[u,b,t]
                set_deactivation(Q_HC[u,b,t] == 0.0 && Q_HC[u,b,t] ==0.0)
# It should be always true that if the internal Temperature is between the Set point values, then both the cooling flux power and the heating flux power should be turn off
            end
        end
    end
end