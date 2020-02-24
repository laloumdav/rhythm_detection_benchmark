#################
### empJTK ###
#################

args = commandArgs(trailingOnly=TRUE)
species <- args[1]
tissue <- args[2]

tissue.file <- paste(tissue, ".txt", sep="")

main.dir <- "~/Documents/rhythm_detection_benchmark"

file.dir <- paste(paste(main.dir, "DATA", species, tissue, sep = "/"), "/", sep="")

# check if working directory exist
if (dir.exists(file.dir)==FALSE){
  stop("Directory/file doesn't exist. Maybe: mouse_RNAseq or mouse_microarray ???")
}

empJTK.python.directory <- paste(main.dir, "/scripts/empJTK_master/", sep = "")
  

# check if empJTK python script is present
if (file.exists(paste(empJTK.python.directory, "eJTK-CalcP.py", sep = ""))==FALSE){
  stop("empJTK files are missing. Should be in rhythm_detection_benchmark/scripts/empJTK_master/")
}



### SCRIPT ###
system(paste("echo Running the file :", tissue.file, sep=" "))

file.name <- paste(file.dir, tissue.file, sep = "")
raw.dataset <- read.table(file.name, head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)


waveforms <- "ref_files/waveform_cosine.txt"
period <- "ref_files/period24.txt"
phases.searched <- "ref_files/phases_00-22_by2.txt"
asymmetries.searched <- "ref_files/asymmetries_02-22_by2.txt"
output <- "empJTK"

### Bash script for empJTK ###
system(paste("chmod 755 ", empJTK.python.directory, "eJTK-CalcP.py", sep=""), wait=TRUE)
system(paste("cd ", empJTK.python.directory, "bin/ ; python setup.py build_ext --inplace", sep=""), wait = TRUE)

empJTK_bash.script <- paste("./eJTK-CalcP.py",
                            "-f", file.name,
                            "-w", waveforms,
                            "-p", period,
                            "-s", phases.searched,
                            "-a", asymmetries.searched,
                            "-x", output,
                            sep = " ")

empJTK_bash.script <- paste("cd ", empJTK.python.directory, " ; ", empJTK_bash.script, sep="")
empJTK_bash.script <- paste("echo waveform searched:", " ; ", "head ", empJTK.python.directory, "ref_files/waveform_cosine.txt", " ; ", empJTK_bash.script, sep="")
empJTK_bash.script <- paste("echo period searched:", " ; ", "head ", empJTK.python.directory, "ref_files/period24.txt", " ; ", empJTK_bash.script, sep="")
system(empJTK_bash.script, wait=TRUE)

# Then we get back the "GAMMA" empJTK files:
empJTK.files <- list.files(file.dir)[grep(tissue, list.files(file.dir))]
empJTK.files <- empJTK.files[grep("empJTK", empJTK.files)]
empJTK.gamma.file <- empJTK.files[grep("Gamma", empJTK.files)]
if (length(empJTK.gamma.file)!=1) { print("ERROR: empJTK_Gamma_output_file doesn't exist or is not unique") }
empJTK.gamma.file <- paste(file.dir, empJTK.gamma.file, sep="")
empJTK <- read.table(empJTK.gamma.file, head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# We keep only info needed
empJTK <- empJTK[,c("ID", "P", "empP", "Period", "Phase", "Max_Amp")]
colnames(empJTK) <- c("ID", "raw.pvalue", "default.pvalue", "period", "phase", "amplitude")
# We remove empJTK_output_files
empJTK.files <- paste(file.dir, empJTK.files, sep = "")
empJTK.files <- paste(empJTK.files, collapse = " ")
remove.empJTK.files_bash.script <- paste("rm", empJTK.files, sep = " ")
system(remove.empJTK.files_bash.script, wait = TRUE)


# Write the new one:
write.table(empJTK, paste(file.dir, "empJTK.txt", sep = ""), row.names = F, quote = F, sep = "\t")
