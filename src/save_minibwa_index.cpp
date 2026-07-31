#include "save_minibwa_index.hpp"
#include "sequence_table_reader.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <cstdlib>
#include <filesystem>
#include <fstream>

namespace duckdb {

// Same mkdtemp/$TMPDIR convention as align_bowtie2::MakeTempIndexDir, staging
// a FASTA for mb_idx_build() rather than an index output directory.
std::string MakeMiniBWATempDir(const char *caller) {
	const char *tmp = std::getenv("TMPDIR");
	if (!tmp || !*tmp) {
		tmp = "/tmp";
	}
	std::string tmpl = std::string(tmp) + "/miint-mb-XXXXXX";
	std::vector<char> buf(tmpl.begin(), tmpl.end());
	buf.push_back('\0');
	if (::mkdtemp(buf.data()) == nullptr) {
		throw IOException("%s: failed to create temp dir under '%s' (errno=%d)", caller, std::string(tmp), errno);
	}
	return std::string(buf.data());
}

void WriteMiniBWASubjectsFasta(const std::string &path, const std::vector<miint::AlignmentSubject> &subjects,
                               const std::string &subject_table, const char *caller) {
	std::ofstream out(path, std::ios::binary);
	if (!out) {
		throw IOException("%s: failed to open temp FASTA for writing: %s", caller, path);
	}
	for (const auto &subject : subjects) {
		if (subject.sequence.empty()) {
			throw InvalidInputException("%s: subject '%s' in table '%s' has an empty sequence", caller,
			                            subject.read_id, subject_table);
		}
		out << '>' << subject.read_id << '\n' << subject.sequence << '\n';
	}
	if (!out) {
		throw IOException("%s: write error while staging temp FASTA: %s", caller, path);
	}
}

unique_ptr<FunctionData> SaveMiniBWAIndexTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                             vector<LogicalType> &return_types,
                                                             vector<std::string> &names) {
	auto data = make_uniq<Data>();

	if (input.inputs.size() < 2) {
		throw BinderException("save_minibwa_index requires subject_table and output_path parameters");
	}
	data->subject_table = input.inputs[0].ToString();
	data->output_path = input.inputs[1].ToString();

	auto sa_bit_param = input.named_parameters.find("sa_bit");
	if (sa_bit_param != input.named_parameters.end() && !sa_bit_param->second.IsNull()) {
		data->config.sa_bit = sa_bit_param->second.GetValue<int32_t>();
		if (data->config.sa_bit <= 0) {
			throw InvalidInputException("sa_bit must be > 0");
		}
	}

	auto subject_schema = ValidateSequenceTableSchema(context, data->subject_table, /*allow_bigint=*/true);
	data->subjects = ReadSubjectTable(context, data->subject_table, subject_schema);

	if (data->subjects.empty()) {
		throw BinderException("Subject table '%s' is empty. Cannot save index for empty subject set.",
		                      data->subject_table);
	}

	for (const auto &name : data->names) {
		names.emplace_back(name);
	}
	for (const auto &type : data->types) {
		return_types.emplace_back(type);
	}

	return data;
}

unique_ptr<GlobalTableFunctionState> SaveMiniBWAIndexTableFunction::InitGlobal(ClientContext &context,
                                                                               TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto gstate = make_uniq<GlobalState>();

	std::string temp_dir = MakeMiniBWATempDir("save_minibwa_index");
	std::string fasta_path = temp_dir + "/subjects.fa";

	try {
		WriteMiniBWASubjectsFasta(fasta_path, data.subjects, data.subject_table, "save_minibwa_index");

		int rc = mb_idx_build(fasta_path.c_str(), data.output_path.c_str(), data.config.sa_bit, /*n_thread=*/1,
		                      /*is_meth=*/0, /*seed=*/11);
		if (rc != 0) {
			throw IOException("mb_idx_build failed (rc=%d) for staged FASTA: %s", rc, fasta_path);
		}
	} catch (...) {
		std::error_code ec;
		std::filesystem::remove_all(temp_dir, ec); // best-effort; can't propagate ec from an exception path
		throw;
	}

	std::error_code ec;
	std::filesystem::remove_all(temp_dir, ec);

	return gstate;
}

unique_ptr<LocalTableFunctionState> SaveMiniBWAIndexTableFunction::InitLocal(ExecutionContext &context,
                                                                             TableFunctionInitInput &input,
                                                                             GlobalTableFunctionState *global_state) {
	return make_uniq<LocalState>();
}

void SaveMiniBWAIndexTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &global_state = data_p.global_state->Cast<GlobalState>();

	if (global_state.done) {
		output.SetCardinality(0);
		return;
	}

	output.data[0].SetValue(0, Value::BOOLEAN(true));
	output.data[1].SetValue(0, Value(bind_data.output_path));
	output.data[2].SetValue(0, Value::BIGINT(static_cast<int64_t>(bind_data.subjects.size())));

	output.SetCardinality(1);
	global_state.done = true;
}

TableFunction SaveMiniBWAIndexTableFunction::GetFunction() {
	auto tf = TableFunction("save_minibwa_index", {LogicalType::VARCHAR, LogicalType::VARCHAR}, Execute, Bind,
	                        InitGlobal, InitLocal);

	tf.named_parameters["sa_bit"] = LogicalType::INTEGER;

	return tf;
}

void SaveMiniBWAIndexTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
