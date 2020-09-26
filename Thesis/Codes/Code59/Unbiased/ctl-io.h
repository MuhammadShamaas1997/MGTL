/* THIS FILE WAS AUTOMATICALLY GENERATED.  DO NOT MODIFY! */
/* generated from the file: meep.scm */

#ifndef CTL_IO_H
#define CTL_IO_H

#include <ctl.h>

#define CXX_CTL_IO 1
namespace ctlio {

/******* Type declarations *******/

typedef struct material_type_struct {
enum { MATERIAL_TYPE_SELF, MATERIAL_FUNCTION, PERFECT_METAL, MEDIUM } which_subclass;
union {
struct material_function_struct *material_function_data;
struct perfect_metal_struct *perfect_metal_data;
struct medium_struct *medium_data;
} subclass;
} material_type;
#define MATERIAL_TYPE_ABSTRACT 1

typedef struct susceptibility_struct {
vector3 sigma_offdiag;
vector3 sigma_diag;
enum { SUSCEPTIBILITY_SELF, MULTILEVEL_ATOM, GYROTROPIC_SATURATED_SUSCEPTIBILITY, DRUDE_SUSCEPTIBILITY, LORENTZIAN_SUSCEPTIBILITY } which_subclass;
union {
struct multilevel_atom_struct *multilevel_atom_data;
struct gyrotropic_saturated_susceptibility_struct *gyrotropic_saturated_susceptibility_data;
struct drude_susceptibility_struct *drude_susceptibility_data;
struct lorentzian_susceptibility_struct *lorentzian_susceptibility_data;
} subclass;
} susceptibility;

typedef struct lorentzian_susceptibility_struct {
number frequency;
number gamma;
enum { LORENTZIAN_SUSCEPTIBILITY_SELF, GYROTROPIC_LORENTZIAN_SUSCEPTIBILITY, NOISY_LORENTZIAN_SUSCEPTIBILITY } which_subclass;
union {
struct gyrotropic_lorentzian_susceptibility_struct *gyrotropic_lorentzian_susceptibility_data;
struct noisy_lorentzian_susceptibility_struct *noisy_lorentzian_susceptibility_data;
} subclass;
} lorentzian_susceptibility;

typedef struct drude_susceptibility_struct {
number frequency;
number gamma;
enum { DRUDE_SUSCEPTIBILITY_SELF, GYROTROPIC_DRUDE_SUSCEPTIBILITY, NOISY_DRUDE_SUSCEPTIBILITY } which_subclass;
union {
struct gyrotropic_drude_susceptibility_struct *gyrotropic_drude_susceptibility_data;
struct noisy_drude_susceptibility_struct *noisy_drude_susceptibility_data;
} subclass;
} drude_susceptibility;

typedef struct noisy_lorentzian_susceptibility_struct {
number noise_amp;
} noisy_lorentzian_susceptibility;

typedef struct noisy_drude_susceptibility_struct {
number noise_amp;
} noisy_drude_susceptibility;

typedef struct gyrotropic_lorentzian_susceptibility_struct {
vector3 bias;
} gyrotropic_lorentzian_susceptibility;

typedef struct gyrotropic_drude_susceptibility_struct {
vector3 bias;
} gyrotropic_drude_susceptibility;

typedef struct gyrotropic_saturated_susceptibility_struct {
number frequency;
number gamma;
vector3 bias;
number alpha;
} gyrotropic_saturated_susceptibility;

typedef struct {
int num_items;
susceptibility *items;
} susceptibility_list;

typedef struct medium_struct {
vector3 epsilon_diag;
vector3 epsilon_offdiag;
vector3 mu_diag;
vector3 mu_offdiag;
susceptibility_list E_susceptibilities;
susceptibility_list H_susceptibilities;
vector3 E_chi2_diag;
vector3 E_chi3_diag;
vector3 H_chi2_diag;
vector3 H_chi3_diag;
vector3 D_conductivity_diag;
vector3 B_conductivity_diag;
} medium;

typedef struct perfect_metal_struct {
} perfect_metal;

typedef struct material_function_struct {
function material_func;
} material_function;

typedef struct geometric_object_struct {
material_type material;
vector3 center;
enum { GEOMETRIC_OBJECT_SELF, PRISM, BLOCK, SPHERE, CYLINDER, COMPOUND_GEOMETRIC_OBJECT } which_subclass;
union {
struct prism_struct *prism_data;
struct block_struct *block_data;
struct sphere_struct *sphere_data;
struct cylinder_struct *cylinder_data;
struct compound_geometric_object_struct *compound_geometric_object_data;
} subclass;
} geometric_object;

typedef struct {
int num_items;
geometric_object *items;
} geometric_object_list;

typedef struct compound_geometric_object_struct {
geometric_object_list component_objects;
} compound_geometric_object;

typedef struct cylinder_struct {
vector3 axis;
number radius;
number height;
enum { CYLINDER_SELF, WEDGE, CONE } which_subclass;
union {
struct wedge_struct *wedge_data;
struct cone_struct *cone_data;
} subclass;
} cylinder;

typedef struct cone_struct {
number radius2;
} cone;

typedef struct wedge_struct {
number wedge_angle;
vector3 wedge_start;
vector3 e1;
vector3 e2;
} wedge;

typedef struct sphere_struct {
number radius;
} sphere;

typedef struct block_struct {
vector3 e1;
vector3 e2;
vector3 e3;
vector3 size;
matrix3x3 projection_matrix;
enum { BLOCK_SELF, ELLIPSOID } which_subclass;
union {
struct ellipsoid_struct *ellipsoid_data;
} subclass;
} block;

typedef struct {
int num_items;
vector3 *items;
} vector3_list;

typedef struct {
int num_items;
number *items;
} number_list;

typedef struct prism_struct {
vector3_list vertices;
number height;
vector3 axis;
vector3_list vertices_p;
vector3 centroid;
number_list workspace;
matrix3x3 m_c2p;
matrix3x3 m_p2c;
} prism;

typedef struct ellipsoid_struct {
vector3 inverse_semi_axes;
} ellipsoid;

typedef struct lattice_struct {
vector3 basis1;
vector3 basis2;
vector3 basis3;
vector3 size;
vector3 basis_size;
vector3 b1;
vector3 b2;
vector3 b3;
matrix3x3 basis;
matrix3x3 metric;
} lattice;

typedef struct transition_struct {
integer from_level;
integer to_level;
number transition_rate;
number frequency;
vector3 sigma_diag;
number gamma;
number pumping_rate;
} transition;

typedef struct {
int num_items;
transition *items;
} transition_list;

typedef struct multilevel_atom_struct {
number_list initial_populations;
transition_list transitions;
} multilevel_atom;

typedef struct symmetry_struct {
integer direction;
cnumber phase;
enum { SYMMETRY_SELF, MIRROR_SYM, ROTATE4_SYM, ROTATE2_SYM } which_subclass;
union {
struct mirror_sym_struct *mirror_sym_data;
struct rotate4_sym_struct *rotate4_sym_data;
struct rotate2_sym_struct *rotate2_sym_data;
} subclass;
} symmetry;

typedef struct rotate2_sym_struct {
} rotate2_sym;

typedef struct rotate4_sym_struct {
} rotate4_sym;

typedef struct mirror_sym_struct {
} mirror_sym;

typedef struct pml_struct {
number thickness;
integer direction;
integer side;
number strength;
number R_asymptotic;
number mean_stretch;
function pml_profile;
enum { PML_SELF, ABSORBER } which_subclass;
union {
struct absorber_struct *absorber_data;
} subclass;
} pml;

typedef struct absorber_struct {
} absorber;

typedef struct volume_class_struct {
vector3 center;
vector3 size;
} volume_class;

typedef struct src_time_struct {
boolean is_integratedp;
enum { SRC_TIME_SELF, CUSTOM_SRC, GAUSSIAN_SRC, CONTINUOUS_SRC } which_subclass;
union {
struct custom_src_struct *custom_src_data;
struct gaussian_src_struct *gaussian_src_data;
struct continuous_src_struct *continuous_src_data;
} subclass;
} src_time;

typedef struct continuous_src_struct {
number frequency;
number start_time;
number end_time;
number width;
number cutoff;
SCM swigval;
} continuous_src;

typedef struct gaussian_src_struct {
number frequency;
number width;
number start_time;
number cutoff;
SCM swigval;
} gaussian_src;

typedef struct custom_src_struct {
function src_func;
number start_time;
number end_time;
SCM swigval;
} custom_src;

typedef struct source_struct {
src_time src;
integer component;
vector3 center;
vector3 size;
cnumber amplitude;
SCM amp_func;
enum { SOURCE_SELF, EIGENMODE_SOURCE } which_subclass;
union {
struct eigenmode_source_struct *eigenmode_source_data;
} subclass;
} source;

typedef struct eigenmode_source_struct {
vector3 eig_lattice_size;
vector3 eig_lattice_center;
integer component;
integer direction;
integer eig_band;
vector3 eig_kpoint;
boolean eig_match_freqp;
integer eig_parity;
integer eig_resolution;
number eig_tolerance;
} eigenmode_source;

typedef struct flux_region_struct {
vector3 center;
vector3 size;
integer direction;
cnumber weight;
} flux_region;

typedef struct energy_region_struct {
vector3 center;
vector3 size;
integer direction;
cnumber weight;
} energy_region;

typedef struct force_region_struct {
vector3 center;
vector3 size;
integer direction;
cnumber weight;
} force_region;

typedef struct near2far_region_struct {
vector3 center;
vector3 size;
integer direction;
cnumber weight;
} near2far_region;

typedef struct {
int num_items;
cvector3 *items;
} cvector3_list;

typedef struct {
int num_items;
cnumber *items;
} cnumber_list;

typedef struct {
int num_items;
material_type *items;
} material_type_list;

typedef struct {
int num_items;
pml *items;
} pml_list;

typedef struct {
int num_items;
symmetry *items;
} symmetry_list;

/******* Input variables *******/
extern char* epsilon_input_file;
extern integer dimensions;
extern material_type default_material;
extern vector3 geometry_center;
extern lattice geometry_lattice;
extern geometric_object_list geometry;
extern boolean ensure_periodicity;

/******* Output variables *******/

extern int num_read_input_vars;
extern int num_write_output_vars;

extern SCM read_input_vars(void);
extern SCM write_output_vars(void);
extern SCM destroy_input_vars(void);
extern SCM destroy_output_vars(void);

/******* external-functions *******/

extern matrix3x3 square_basis(matrix3x3, vector3);
extern SCM square_basis_aux(SCM arg_scm_0, SCM arg_scm_1);

extern number range_overlap_with_object(vector3, vector3, geometric_object, number, integer);
extern SCM range_overlap_with_object_aux(SCM arg_scm_0, SCM arg_scm_1, SCM arg_scm_2, SCM arg_scm_3, SCM arg_scm_4);

extern void display_geometric_object_info(integer, geometric_object);
extern SCM display_geometric_object_info_aux(SCM arg_scm_0, SCM arg_scm_1);

extern boolean point_in_periodic_objectp(vector3, geometric_object);
extern SCM point_in_periodic_objectp_aux(SCM arg_scm_0, SCM arg_scm_1);

extern vector3 normal_to_object(vector3, geometric_object);
extern SCM normal_to_object_aux(SCM arg_scm_0, SCM arg_scm_1);

extern boolean point_in_objectp(vector3, geometric_object);
extern SCM point_in_objectp_aux(SCM arg_scm_0, SCM arg_scm_1);


extern void export_external_functions(void);

/******* class input function prototypes *******/

extern void near2far_region_input(SCM so, near2far_region *o);
extern void force_region_input(SCM so, force_region *o);
extern void energy_region_input(SCM so, energy_region *o);
extern void flux_region_input(SCM so, flux_region *o);
extern void eigenmode_source_input(SCM so, eigenmode_source *o);
extern void source_input(SCM so, source *o);
extern void custom_src_input(SCM so, custom_src *o);
extern void gaussian_src_input(SCM so, gaussian_src *o);
extern void continuous_src_input(SCM so, continuous_src *o);
extern void src_time_input(SCM so, src_time *o);
extern void volume_class_input(SCM so, volume_class *o);
extern void absorber_input(SCM so, absorber *o);
extern void pml_input(SCM so, pml *o);
extern void mirror_sym_input(SCM so, mirror_sym *o);
extern void rotate4_sym_input(SCM so, rotate4_sym *o);
extern void rotate2_sym_input(SCM so, rotate2_sym *o);
extern void symmetry_input(SCM so, symmetry *o);
extern void multilevel_atom_input(SCM so, multilevel_atom *o);
extern void transition_input(SCM so, transition *o);
extern void lattice_input(SCM so, lattice *o);
extern void ellipsoid_input(SCM so, ellipsoid *o);
extern void prism_input(SCM so, prism *o);
extern void block_input(SCM so, block *o);
extern void sphere_input(SCM so, sphere *o);
extern void wedge_input(SCM so, wedge *o);
extern void cone_input(SCM so, cone *o);
extern void cylinder_input(SCM so, cylinder *o);
extern void compound_geometric_object_input(SCM so, compound_geometric_object *o);
extern void geometric_object_input(SCM so, geometric_object *o);
extern void material_function_input(SCM so, material_function *o);
extern void perfect_metal_input(SCM so, perfect_metal *o);
extern void medium_input(SCM so, medium *o);
extern void gyrotropic_saturated_susceptibility_input(SCM so, gyrotropic_saturated_susceptibility *o);
extern void gyrotropic_drude_susceptibility_input(SCM so, gyrotropic_drude_susceptibility *o);
extern void gyrotropic_lorentzian_susceptibility_input(SCM so, gyrotropic_lorentzian_susceptibility *o);
extern void noisy_drude_susceptibility_input(SCM so, noisy_drude_susceptibility *o);
extern void noisy_lorentzian_susceptibility_input(SCM so, noisy_lorentzian_susceptibility *o);
extern void drude_susceptibility_input(SCM so, drude_susceptibility *o);
extern void lorentzian_susceptibility_input(SCM so, lorentzian_susceptibility *o);
extern void susceptibility_input(SCM so, susceptibility *o);
extern void material_type_input(SCM so, material_type *o);

/******* class copy function prototypes *******/

extern void near2far_region_copy(const near2far_region *o0,near2far_region *o);
extern void force_region_copy(const force_region *o0,force_region *o);
extern void energy_region_copy(const energy_region *o0,energy_region *o);
extern void flux_region_copy(const flux_region *o0,flux_region *o);
extern void eigenmode_source_copy(const eigenmode_source *o0,eigenmode_source *o);
extern void source_copy(const source *o0,source *o);
extern void custom_src_copy(const custom_src *o0,custom_src *o);
extern void gaussian_src_copy(const gaussian_src *o0,gaussian_src *o);
extern void continuous_src_copy(const continuous_src *o0,continuous_src *o);
extern void src_time_copy(const src_time *o0,src_time *o);
extern void volume_class_copy(const volume_class *o0,volume_class *o);
extern void absorber_copy(const absorber *o0,absorber *o);
extern void pml_copy(const pml *o0,pml *o);
extern void mirror_sym_copy(const mirror_sym *o0,mirror_sym *o);
extern void rotate4_sym_copy(const rotate4_sym *o0,rotate4_sym *o);
extern void rotate2_sym_copy(const rotate2_sym *o0,rotate2_sym *o);
extern void symmetry_copy(const symmetry *o0,symmetry *o);
extern void multilevel_atom_copy(const multilevel_atom *o0,multilevel_atom *o);
extern void transition_copy(const transition *o0,transition *o);
extern void lattice_copy(const lattice *o0,lattice *o);
extern void ellipsoid_copy(const ellipsoid *o0,ellipsoid *o);
extern void prism_copy(const prism *o0,prism *o);
extern void block_copy(const block *o0,block *o);
extern void sphere_copy(const sphere *o0,sphere *o);
extern void wedge_copy(const wedge *o0,wedge *o);
extern void cone_copy(const cone *o0,cone *o);
extern void cylinder_copy(const cylinder *o0,cylinder *o);
extern void compound_geometric_object_copy(const compound_geometric_object *o0,compound_geometric_object *o);
extern void geometric_object_copy(const geometric_object *o0,geometric_object *o);
extern void material_function_copy(const material_function *o0,material_function *o);
extern void perfect_metal_copy(const perfect_metal *o0,perfect_metal *o);
extern void medium_copy(const medium *o0,medium *o);
extern void gyrotropic_saturated_susceptibility_copy(const gyrotropic_saturated_susceptibility *o0,gyrotropic_saturated_susceptibility *o);
extern void gyrotropic_drude_susceptibility_copy(const gyrotropic_drude_susceptibility *o0,gyrotropic_drude_susceptibility *o);
extern void gyrotropic_lorentzian_susceptibility_copy(const gyrotropic_lorentzian_susceptibility *o0,gyrotropic_lorentzian_susceptibility *o);
extern void noisy_drude_susceptibility_copy(const noisy_drude_susceptibility *o0,noisy_drude_susceptibility *o);
extern void noisy_lorentzian_susceptibility_copy(const noisy_lorentzian_susceptibility *o0,noisy_lorentzian_susceptibility *o);
extern void drude_susceptibility_copy(const drude_susceptibility *o0,drude_susceptibility *o);
extern void lorentzian_susceptibility_copy(const lorentzian_susceptibility *o0,lorentzian_susceptibility *o);
extern void susceptibility_copy(const susceptibility *o0,susceptibility *o);
extern void material_type_copy(const material_type *o0,material_type *o);

/******* class equal function prototypes *******/

extern boolean near2far_region_equal(const near2far_region *o0, const near2far_region *o);
extern boolean force_region_equal(const force_region *o0, const force_region *o);
extern boolean energy_region_equal(const energy_region *o0, const energy_region *o);
extern boolean flux_region_equal(const flux_region *o0, const flux_region *o);
extern boolean eigenmode_source_equal(const eigenmode_source *o0, const eigenmode_source *o);
extern boolean source_equal(const source *o0, const source *o);
extern boolean custom_src_equal(const custom_src *o0, const custom_src *o);
extern boolean gaussian_src_equal(const gaussian_src *o0, const gaussian_src *o);
extern boolean continuous_src_equal(const continuous_src *o0, const continuous_src *o);
extern boolean src_time_equal(const src_time *o0, const src_time *o);
extern boolean volume_class_equal(const volume_class *o0, const volume_class *o);
extern boolean absorber_equal(const absorber *o0, const absorber *o);
extern boolean pml_equal(const pml *o0, const pml *o);
extern boolean mirror_sym_equal(const mirror_sym *o0, const mirror_sym *o);
extern boolean rotate4_sym_equal(const rotate4_sym *o0, const rotate4_sym *o);
extern boolean rotate2_sym_equal(const rotate2_sym *o0, const rotate2_sym *o);
extern boolean symmetry_equal(const symmetry *o0, const symmetry *o);
extern boolean multilevel_atom_equal(const multilevel_atom *o0, const multilevel_atom *o);
extern boolean transition_equal(const transition *o0, const transition *o);
extern boolean lattice_equal(const lattice *o0, const lattice *o);
extern boolean ellipsoid_equal(const ellipsoid *o0, const ellipsoid *o);
extern boolean prism_equal(const prism *o0, const prism *o);
extern boolean block_equal(const block *o0, const block *o);
extern boolean sphere_equal(const sphere *o0, const sphere *o);
extern boolean wedge_equal(const wedge *o0, const wedge *o);
extern boolean cone_equal(const cone *o0, const cone *o);
extern boolean cylinder_equal(const cylinder *o0, const cylinder *o);
extern boolean compound_geometric_object_equal(const compound_geometric_object *o0, const compound_geometric_object *o);
extern boolean geometric_object_equal(const geometric_object *o0, const geometric_object *o);
extern boolean material_function_equal(const material_function *o0, const material_function *o);
extern boolean perfect_metal_equal(const perfect_metal *o0, const perfect_metal *o);
extern boolean medium_equal(const medium *o0, const medium *o);
extern boolean gyrotropic_saturated_susceptibility_equal(const gyrotropic_saturated_susceptibility *o0, const gyrotropic_saturated_susceptibility *o);
extern boolean gyrotropic_drude_susceptibility_equal(const gyrotropic_drude_susceptibility *o0, const gyrotropic_drude_susceptibility *o);
extern boolean gyrotropic_lorentzian_susceptibility_equal(const gyrotropic_lorentzian_susceptibility *o0, const gyrotropic_lorentzian_susceptibility *o);
extern boolean noisy_drude_susceptibility_equal(const noisy_drude_susceptibility *o0, const noisy_drude_susceptibility *o);
extern boolean noisy_lorentzian_susceptibility_equal(const noisy_lorentzian_susceptibility *o0, const noisy_lorentzian_susceptibility *o);
extern boolean drude_susceptibility_equal(const drude_susceptibility *o0, const drude_susceptibility *o);
extern boolean lorentzian_susceptibility_equal(const lorentzian_susceptibility *o0, const lorentzian_susceptibility *o);
extern boolean susceptibility_equal(const susceptibility *o0, const susceptibility *o);
extern boolean material_type_equal(const material_type *o0, const material_type *o);

/******* class destruction function prototypes *******/

extern void near2far_region_destroy(near2far_region o);
extern void force_region_destroy(force_region o);
extern void energy_region_destroy(energy_region o);
extern void flux_region_destroy(flux_region o);
extern void eigenmode_source_destroy(eigenmode_source o);
extern void source_destroy(source o);
extern void custom_src_destroy(custom_src o);
extern void gaussian_src_destroy(gaussian_src o);
extern void continuous_src_destroy(continuous_src o);
extern void src_time_destroy(src_time o);
extern void volume_class_destroy(volume_class o);
extern void absorber_destroy(absorber o);
extern void pml_destroy(pml o);
extern void mirror_sym_destroy(mirror_sym o);
extern void rotate4_sym_destroy(rotate4_sym o);
extern void rotate2_sym_destroy(rotate2_sym o);
extern void symmetry_destroy(symmetry o);
extern void multilevel_atom_destroy(multilevel_atom o);
extern void transition_destroy(transition o);
extern void lattice_destroy(lattice o);
extern void ellipsoid_destroy(ellipsoid o);
extern void prism_destroy(prism o);
extern void block_destroy(block o);
extern void sphere_destroy(sphere o);
extern void wedge_destroy(wedge o);
extern void cone_destroy(cone o);
extern void cylinder_destroy(cylinder o);
extern void compound_geometric_object_destroy(compound_geometric_object o);
extern void geometric_object_destroy(geometric_object o);
extern void material_function_destroy(material_function o);
extern void perfect_metal_destroy(perfect_metal o);
extern void medium_destroy(medium o);
extern void gyrotropic_saturated_susceptibility_destroy(gyrotropic_saturated_susceptibility o);
extern void gyrotropic_drude_susceptibility_destroy(gyrotropic_drude_susceptibility o);
extern void gyrotropic_lorentzian_susceptibility_destroy(gyrotropic_lorentzian_susceptibility o);
extern void noisy_drude_susceptibility_destroy(noisy_drude_susceptibility o);
extern void noisy_lorentzian_susceptibility_destroy(noisy_lorentzian_susceptibility o);
extern void drude_susceptibility_destroy(drude_susceptibility o);
extern void lorentzian_susceptibility_destroy(lorentzian_susceptibility o);
extern void susceptibility_destroy(susceptibility o);
extern void material_type_destroy(material_type o);


} /* namespace */

#endif                          /* CTL_IO_H */

