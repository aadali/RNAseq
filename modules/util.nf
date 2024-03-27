def get_sample_number(sample_info){
    def a = file(sample_info, checkIfExists: true).readLines().findAll{it -> !it.startsWith("#")} // the header line startsWith("sample") is true
    return a.size() - 1
}