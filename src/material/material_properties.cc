#include "material_properties.h"

MaterialProperties::MaterialProperties(const std::unordered_map<int, Material>& materials,
                                        bool is_eigen_problem,
                                        bool do_nda,
                                        int number_of_groups,
                                        int number_of_materials,
                                        const std::unordered_set<int>& fissile_ids /* = {} */)
    : materials_(materials),
      is_eigen_problem_(is_eigen_problem),
      do_nda_(do_nda),
      n_group_(number_of_groups),
      n_material_(number_of_materials),
      fissile_ids_(fissile_ids) {
    CheckFissileIDs();
    PopulateFissileMap(); // generates is_material_fissile_, which is used by CheckValidEach and PopulateData

    CheckNumberOfMaterials();
    CheckValidEach();
    CheckNumberOfGroups();
    CheckConsistent();

    PopulateData();
}

MaterialProperties::MaterialProperties(dealii::ParameterHandler& prm)
    : MaterialProperties(
        ParseMaterials(ReadMaterialFileNames(prm)),
        prm.get_bool("do eigenvalue calculations"),
        prm.get_bool("do NDA"),
        prm.get_integer("number of groups"),
        prm.get_integer("number of materials"),
        ReadFissileIDs(prm)) {}

MaterialProperties::~MaterialProperties() {}

std::unordered_map<int, std::string> MaterialProperties::ReadMaterialFileNames(dealii::ParameterHandler& prm) {
  std::unordered_map<int, std::string> result;
  prm.enter_subsection("material ID map");
  const std::vector<std::string> pair_strings = dealii::Utilities::split_string_list(prm.get("material id file name map"), ",");
  for (const std::string& pair_string : pair_strings) {
    const std::vector<std::string> split_pair = dealii::Utilities::split_string_list(pair_string, ':');
    result[std::stoi(split_pair[0])] = split_pair[1];
  }
  prm.leave_subsection();
  return result;
}

std::unordered_set<int> MaterialProperties::ReadFissileIDs(dealii::ParameterHandler& prm) {
  std::unordered_set<int> result;
  prm.enter_subsection("fissile material IDs");
  const std::vector<std::string> id_strings = dealii::Utilities::split_string_list(prm.get("fissile material ids"), ",");
  for (const std::string& id_string : id_strings) {
    result.insert(std::stoi(id_string));
  }
  prm.leave_subsection();
  return result;
}

std::unordered_map<int, Material> MaterialProperties::ParseMaterials(const std::unordered_map<int, std::string>& file_name_map) {
  /*
    Tries serialized format first and then readable text format if parsing as serialized format fails.
    This is because by default, parsing as text format will print errors if it doesn't work, but
    parsing as serialized format will not.
  */
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  std::unordered_map<int, Material> result;
  for (const std::pair<int, std::string>& id_file_pair : file_name_map) {
    result[id_file_pair.first] = Material();
    bool text_format_success = false;
    bool serialized_success = false;

    std::ifstream serialized_ifstream(id_file_pair.second, std::ios::binary);
    AssertThrow(serialized_ifstream.is_open(),
    FailedToFindMaterialFile(id_file_pair.second, id_file_pair.first));
    serialized_success = result[id_file_pair.first].ParseFromIstream(&serialized_ifstream);
    serialized_ifstream.close();

    if (!serialized_success) {
      std::ifstream text_format_ifstream(id_file_pair.second);
      AssertThrow(text_format_ifstream.is_open(),
      FailedToFindMaterialFile(id_file_pair.second, id_file_pair.first));
      google::protobuf::io::IstreamInputStream zero_copy_stream(&text_format_ifstream); //TextFormat::Parse requires different stream type
      text_format_success = google::protobuf::TextFormat::Parse(&zero_copy_stream, &result[id_file_pair.first]);
      text_format_ifstream.close();
    }

    AssertThrow(text_format_success || serialized_success,
      FailedToParseMaterialFile(id_file_pair.second, id_file_pair.first));
  }
  return result;
}

void MaterialProperties::PopulateFissileMap() {
  for (const std::pair<int, Material>& mat_pair : materials_) {
    const int& id = mat_pair.first;
    is_material_fissile_[id] = (fissile_ids_.count(id) > 0);
  }
}

void MaterialProperties::PopulateData() {
  /* 
    unordered maps for fixed source data will be left empty for an eigen problem and
    materials not labeled fissile will be excluded from the unordered_maps for fission data
    Q is assumed to be zero for all groups if not given in the Material
  */

  for (const std::pair<int, Material>& mat_pair : materials_) {
    const int& id = mat_pair.first;
    const Material& material = mat_pair.second;
    const auto vector_props = GetVectorProperties(material);

    sigt_[id] = vector_props.at(Material::SIGMA_T);
    sigs_[id] = GetScatteringMatrix(material);

    if (is_material_fissile_[id]) {
      chi_[id] = vector_props.at(Material::CHI);
      nusigf_[id] = vector_props.at(Material::NU_SIG_F);
    }

    if (!is_eigen_problem_) {
      if (vector_props.count(Material::Q) > 0) {
        q_[id] = vector_props.at(Material::Q);
      }
      else {
        q_[id] = std::vector<double>(n_group_, 0);
      }
    }
  }

  for (const int& id : fissile_ids_) {
    dealii::FullMatrix nusigf_column_matrix(n_group_, 1, nusigf_.at(id).data());
    dealii::FullMatrix chi_row_matrix(1, n_group_, chi_.at(id).data());
    chi_nusigf_[id] = dealii::FullMatrix<double>(n_group_, n_group_);
    nusigf_column_matrix.mmult(chi_nusigf_[id], chi_row_matrix); // chi_nusigf_[id] = nusigf_column_matrix * chi_row_matrix
  }

  for (const std::pair<int, std::vector<double>>& sigt_pair : sigt_) {
    inv_sigt_[sigt_pair.first] = sigt_pair.second;
    for (double& val : inv_sigt_.at(sigt_pair.first)) {
      val = 1 / val;
    }
  }

  for (const std::pair<int, std::vector<double>>& q_pair : q_) {
    q_per_ster_[q_pair.first] = q_pair.second;
    for (double& val : q_per_ster_.at(q_pair.first)) {
      val *= bconst::kInvFourPi;
    }
  }

  for (const std::pair<int, dealii::FullMatrix<double>>& sigs_pair : sigs_) {
    sigs_per_ster_[sigs_pair.first] = sigs_pair.second;
    sigs_per_ster_[sigs_pair.first] *= bconst::kInvFourPi;
  }

  for (const std::pair<int, dealii::FullMatrix<double>>& chi_nusigf_pair : chi_nusigf_) {
    chi_nusigf_per_ster_[chi_nusigf_pair.first] = chi_nusigf_pair.second;
    chi_nusigf_per_ster_[chi_nusigf_pair.first] *= bconst::kInvFourPi;
  }
}

void MaterialProperties::CheckFissileIDs() const {
  AssertThrow(!(is_eigen_problem_ && fissile_ids_.empty()), NoFissileIDs());

  for (const int& id : fissile_ids_) {
    AssertThrow(materials_.count(id) > 0, FissileIDInvalid(id));
  }
}

void MaterialProperties::CheckNumberOfMaterials() const {
  AssertThrow(n_material_ >= 0 && (unsigned int)(n_material_) == materials_.size(),
    WrongNumberOfMaterials(materials_.size(), n_material_));
}

void MaterialProperties::CheckNumberOfGroups() const {
  for (const std::pair<int, Material>& mat_pair : materials_) {
    const Material& mat = mat_pair.second;
    AssertThrow(n_group_ >= 0 && (unsigned int)(n_group_) == mat.number_of_groups(),
      WrongNumberOfGroups(CombinedName(mat_pair), mat.number_of_groups(), n_group_));
  }
}

void MaterialProperties::CheckValid(const Material& material, const bool require_fission_data /* = false*/, std::string name /* = ""*/) {
  /* 
    calls GetVectorProperties and GetScatteringMatrix, which can throw exceptions, 
    but the order of checks is such that they shouldn't
  */
  if (name == "") {
    name = CombinedName(material);
  }

  // MultipleDefinition
  std::unordered_set<Material::VectorId> vector_ids_in_material;

  for (const Material_VectorProperty& vec_prop : material.vector_property()) {
    AssertThrow(vector_ids_in_material.count(vec_prop.id()) == 0 || vec_prop.id() == Material::UNKNOWN_VECTOR,
      MultipleDefinition(name, Material_VectorId_descriptor()->value(vec_prop.id())->name()));

    vector_ids_in_material.insert(vec_prop.id());
  }

  std::unordered_set<Material::MatrixId> matrix_ids_in_material;

  for (const Material_MatrixProperty& mat_prop : material.matrix_property()) {
    AssertThrow(matrix_ids_in_material.count(mat_prop.id()) == 0 || mat_prop.id() == Material::UNKNOWN_MATRIX,
      MultipleDefinition(name, Material_MatrixId_descriptor()->value(mat_prop.id())->name()));

    matrix_ids_in_material.insert(mat_prop.id());
  }

  // MissingProperty
  std::unordered_set<Material::VectorId> required_vector_props = {Material::ENERGY_GROUPS, Material::SIGMA_T};

  if (require_fission_data) {
    required_vector_props.insert(Material::NU_SIG_F);
    required_vector_props.insert(Material::CHI);
  }

  std::unordered_map<Material::VectorId, std::vector<double>> vector_props = GetVectorProperties(material);

  for (Material::VectorId id : required_vector_props) {
    AssertThrow(vector_props.count(id) == 1,
      MissingProperty(name, Material_VectorId_descriptor()->value(id)->name()));
  }

  AssertThrow(matrix_ids_in_material.count(Material::SIGMA_S) > 0,
    MissingProperty(name, Material_MatrixId_descriptor()->value(Material::SIGMA_S)->name()));

  // WrongNumberOfValues
  if (vector_props.count(Material::Q) > 0) {
    required_vector_props.insert(Material::Q);
  }
  const unsigned int& n = material.number_of_groups();
  std::unordered_map<Material::VectorId, unsigned int> required_count =
   {{Material::ENERGY_GROUPS, n+1}, {Material::SIGMA_T, n}, {Material::Q, n}, {Material::NU_SIG_F, n}, {Material::CHI, n}};

  for (Material::VectorId id : required_vector_props) {
    AssertThrow(vector_props.at(id).size() == required_count.at(id),
      WrongNumberOfValues(name, Material_VectorId_descriptor()->value(id)->name(),
        vector_props.at(id).size(), n));
  }

  for (const Material_MatrixProperty& mat_prop : material.matrix_property()) {
    if (mat_prop.id() == Material::SIGMA_S) {
      AssertThrow((unsigned int)mat_prop.value().size() == material.number_of_groups()*material.number_of_groups(),
          WrongNumberOfValues(name, Material_MatrixId_descriptor()->value(Material::SIGMA_S)->name(),
            mat_prop.value().size(), material.number_of_groups()));
    }
  }

  // WrongSign
  const std::unordered_set<Material::VectorId> required_non_negative = required_vector_props;
  for (Material::VectorId id : required_non_negative) {
    for (unsigned int i = 0; i < vector_props.at(id).size(); ++i) {
      AssertThrow(vector_props.at(id)[i] >= 0,
        WrongSign(name, Material_VectorId_descriptor()->value(id)->name(),
          vector_props.at(id)[i], i));
    }
  }

  const dealii::FullMatrix<double> scattering_matrix = GetScatteringMatrix(material);
  for (auto iter = scattering_matrix.begin(); iter < scattering_matrix.end(); ++iter) {
    AssertThrow(iter->value() >= 0,
      WrongSign(name, Material_MatrixId_descriptor()->value(Material::SIGMA_S)->name(),
        iter->value(), iter->row()*n + iter->column()));
  }

  // EnergyGroupBoundariesNotSorted
  for (auto i = vector_props.at(Material::ENERGY_GROUPS).begin()+1; i < vector_props.at(Material::ENERGY_GROUPS).end()-1; ++i) {
    // require that every set of three consecutive values is in either strictly ascending or strictly descending order
    AssertThrow((*(i-1) > *i && *i > *(i+1)) || (*(i-1) < *i && *i < *(i+1)),
      EnergyGroupBoundariesNotSorted(name, *(i-1), *i, *(i+1)));
  }

  // ChiDoesNotSumToOne
  if (required_vector_props.count(Material::CHI) > 0) {
    // allow a maximum error from 1 of two units in the last place
    AssertThrow(abs(PreciseSum(vector_props.at(Material::CHI)) - 1) <= 2*std::numeric_limits<double>::epsilon(),
      ChiDoesNotSumToOne(name, PreciseSum(vector_props.at(Material::CHI))));
  }
}

void MaterialProperties::CheckValidEach() const {
  for (const std::pair<int, Material>& mat_pair : materials_) {
    CheckValid(mat_pair.second, is_material_fissile_.at(mat_pair.first), CombinedName(mat_pair));
  }
}

void MaterialProperties::CheckConsistent(const std::unordered_map<int, Material>& materials) {
  if (materials.empty()) {
    return;
  }

  const std::vector<double> first_mat_energies = GetVectorProperty(materials.cbegin()->second, Material::ENERGY_GROUPS);

  for (std::unordered_map<int, Material>::const_iterator iter = materials.cbegin(); iter != materials.cend(); ++iter) {
    AssertThrow(materials.cbegin()->second.number_of_groups() == iter->second.number_of_groups(),
      NumberOfGroupsMismatch(CombinedName(*materials.cbegin()), CombinedName(*iter)));

    const std::vector<double> this_mat_energies = GetVectorProperty(iter->second, Material::ENERGY_GROUPS);

    AssertThrow(first_mat_energies.size() == this_mat_energies.size(),
      EnergyGroupsMismatch(CombinedName(*materials.cbegin()), CombinedName(*iter),
        std::min(first_mat_energies.size(), this_mat_energies.size())+1));

    for (unsigned int i = 0; i < first_mat_energies.size(); ++i) {
      AssertThrow(first_mat_energies[i] == this_mat_energies[i],
        EnergyGroupsMismatch(CombinedName(*materials.cbegin()), CombinedName(*iter), i));
    }
  }
}

void MaterialProperties::CheckConsistent() const {
  CheckConsistent(materials_);
}

std::string MaterialProperties::CombinedName(const std::pair<int, Material>& id_material_pair, const std::string& delimiter /* = "|" */) {
  const std::string q = "\"";
  const std::string& d = delimiter;
  const Material& mat = id_material_pair.second;
  std::string format = "number" + d + "full_name" + d + "abbreviation" + d + "id";
  std::string info = std::to_string(id_material_pair.first) + d + q+mat.full_name()+q + d + q+mat.abbreviation()+q + d + q+mat.id()+q;
  return format + " = " + info;
}

std::string MaterialProperties::CombinedName(const Material& material, const std::string& delimiter /* = "|" */) {
  const std::string q = "\"";
  const std::string& d = delimiter;
  std::string format = "full_name" + d + "abbreviation" + d + "id";
  std::string info = q+material.full_name()+q + d + q+material.abbreviation()+q + d + q+material.id()+q;
  return format + " = " + info;
}

std::unordered_map<Material::VectorId, std::vector<double>> MaterialProperties::GetVectorProperties(const Material& material) {
  std::unordered_map<Material::VectorId, std::vector<double>> result;
  for (const Material_VectorProperty& vec_prop : material.vector_property()) {
    AssertThrow(result.count(vec_prop.id()) == 0 || vec_prop.id() == Material::UNKNOWN_VECTOR,
      MultipleDefinition(CombinedName(material), Material_VectorId_descriptor()->value(vec_prop.id())->name()));

    std::vector<double> new_vector(vec_prop.value().cbegin(), vec_prop.value().cend());
    result[vec_prop.id()] = new_vector;
  }
  return result;
}

std::vector<double> MaterialProperties::GetVectorProperty(const Material& material, Material::VectorId property_id) {
  for (const Material_VectorProperty& vec_prop : material.vector_property()) {
    if (vec_prop.id() == property_id) {
      return std::vector<double>(vec_prop.value().cbegin(), vec_prop.value().cend());
    }
  }
  return {};
}

dealii::FullMatrix<double> MaterialProperties::GetScatteringMatrix(const Material& material) {
  for (const Material_MatrixProperty& mat_prop : material.matrix_property()) {
    if (mat_prop.id() == Material::SIGMA_S) {
      std::vector<double> raw_values(mat_prop.value().cbegin(), mat_prop.value().cend());
      
      AssertThrow(raw_values.size() == material.number_of_groups()*material.number_of_groups(),
        WrongNumberOfValues(CombinedName(material), Material_MatrixId_descriptor()->value(Material::SIGMA_S)->name(),
          raw_values.size(), material.number_of_groups()));

      dealii::FullMatrix<double> scattering_matrix(material.number_of_groups(), material.number_of_groups());
      scattering_matrix.fill(raw_values.data());
      return scattering_matrix;
    }
  }
  return dealii::FullMatrix<double>();
}

double MaterialProperties::PreciseSum(const std::vector<double>& values) {
  // uses Kahan summation algorithm
  double sum = 0;
  double error = 0;
  for (const double& addend : values) {
    const double corrected_addend = addend - error;
    const double new_sum = sum + corrected_addend;
    error = (new_sum - sum) - corrected_addend; // holds low order digits of corrected_addend lost when making new_sum
    sum = new_sum;
  }
  return sum;
}

std::unordered_map<int, bool> MaterialProperties::GetFissileIDMap() const {
  return is_material_fissile_;
}

std::unordered_map<int, std::vector<double>> MaterialProperties::GetSigT() const {
  return sigt_;
}

std::unordered_map<int, std::vector<double>> MaterialProperties::GetInvSigT() const {
  return inv_sigt_;
}

std::unordered_map<int, std::vector<double>> MaterialProperties::GetQ() const {
  return q_;
}

std::unordered_map<int, std::vector<double>> MaterialProperties::GetQPerSter() const {
  return q_per_ster_;
}

std::unordered_map<int, std::vector<double>> MaterialProperties::GetNuSigF() const {
  return nusigf_;
}

std::unordered_map<int, dealii::FullMatrix<double>> MaterialProperties::GetSigS() const {
  return sigs_;
}

std::unordered_map<int, dealii::FullMatrix<double>> MaterialProperties::GetSigSPerSter() const {
  return sigs_per_ster_;
}

std::unordered_map<int, dealii::FullMatrix<double>> MaterialProperties::GetChiNuSigF() const {
  return chi_nusigf_;
}

std::unordered_map<int, dealii::FullMatrix<double>> MaterialProperties::GetChiNuSigFPerSter() const {
  return chi_nusigf_per_ster_;
}
