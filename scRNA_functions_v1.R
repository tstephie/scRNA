# scRNA Functions v1
# 06/17/19

create_object <- function(object_name, directory, project_name, min.cells, min.features) {
  assign(object_name, CreateSeuratObject(counts = Read10X(directory), project = project_name, 
                                         min.cells = min.cells, min.features = min.features))
}