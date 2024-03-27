
max_cpus = "nproc".execute().text.toInteger()
println(Math.floor(max_cpus/5).toInteger())
println(max_cpus / 5)