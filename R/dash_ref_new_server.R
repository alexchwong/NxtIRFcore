server_ref_new <- function(id, refresh_tab, volumes) {
    moduleServer(id, function(input, output, session) {
        settings_newref <- setreactive_newref()

        observe({  
            shinyDirChoose(input, "dir_reference_path", 
                roots = volumes(), session = session)
            output$txt_reference_path <- renderText({
                validate(need(input$dir_reference_path, 
                    "Please select reference path"))
                settings_newref$newref_path = 
                    parseDirPath(volumes(), input$dir_reference_path)
            })
        })
        observe({
            shinyFileChoose(input, "file_genome", 
                roots = volumes(), 
                session = session, filetypes = c("fa", "fasta", "gz"))
            if(!is.null(input$file_genome)){
                file_selected<-parseFilePaths(volumes(), input$file_genome)
                settings_newref$newref_fasta = 
                    as.character(file_selected$datapath)
                output$txt_genome <- renderText(
                    as.character(file_selected$datapath))
            }
        })
        observe({  
            shinyFileChoose(input, "file_gtf", roots = volumes(), 
                session = session, filetypes = c("gtf", "gz"))
            if(!is.null(input$file_gtf)){
                file_selected<-parseFilePaths(volumes(), input$file_gtf)
                settings_newref$newref_gtf = 
                    as.character(file_selected$datapath)
                output$txt_gtf <- renderText(
                    as.character(file_selected$datapath))
            }
        })
        observe({  
            shinyFileChoose(input, "file_mappa", roots = volumes(), 
                session = session, filetypes = c("txt", "gz"))
            if(!is.null(input$file_mappa)){
                    updateSelectInput(session = session, 
                        inputId = "newref_genome_type", 
                        choices = c("(custom)", "hg38", "mm10", "hg19", "mm9"),
                        selected = "(custom)"
                    )
                file_selected<-parseFilePaths(volumes(), input$file_mappa)
                settings_newref$newref_mappa = 
                    as.character(file_selected$datapath)
            }
        })
        observeEvent(settings_newref$newref_mappa, {
            output$txt_mappa <- renderText(settings_newref$newref_mappa)
        })
        observeEvent(input$clear_mappa, {
            req(input$clear_mappa)
            settings_newref$newref_mappa = ""
        })    
        observe({  
            shinyFileChoose(input, "file_NPA", roots = volumes(), 
                session = session, filetypes = c("bed", "txt", "gz"))
            if(!is.null(input$file_NPA)){
                updateSelectInput(session = session, 
                    inputId = "newref_genome_type", 
                    choices = c("(custom)", "hg38", "mm10", "hg19", "mm9"),
                    selected = "(custom)"
                )
                file_selected<-parseFilePaths(volumes(), input$file_NPA)
                settings_newref$newref_NPA = 
                    as.character(file_selected$datapath)
            }
        })
        observeEvent(settings_newref$newref_NPA, {
            output$txt_NPA <- renderText(settings_newref$newref_NPA)    
        })
        observeEvent(input$clear_NPA, {
            req(input$clear_NPA)
            settings_newref$newref_NPA = ""
        })    
        observe({  
            shinyFileChoose(input, "file_bl", roots = volumes(), 
                session = session, filetypes = c("bed", "txt", "gz"))
            if(!is.null(input$file_bl)){
                file_selected<-parseFilePaths(volumes(), input$file_bl)
                settings_newref$newref_bl = as.character(file_selected$datapath)
            }
        })
        observeEvent(settings_newref$newref_bl, {
            output$txt_bl <- renderText(settings_newref$newref_bl)
        })
        observeEvent(input$clear_bl, {
        req(input$clear_bl)
            settings_newref$newref_bl = ""
        })
        observeEvent(input$newref_genome_type, {
            req(input$newref_genome_type)

            if(input$newref_genome_type == "hg38") {
                settings_newref$newref_NPA = system.file(
                    "extra-input-files/Human_hg38_nonPolyA_ROI.bed", 
                    package = "NxtIRF")
                settings_newref$newref_mappa = 
                    GetMappabilityRef("hg38")
            
            } else if(input$newref_genome_type == "hg19")  {
                settings_newref$newref_NPA = system.file(
                    "extra-input-files/Human_hg19_nonPolyA_ROI.bed", 
                    package = "NxtIRF")
                settings_newref$newref_mappa = 
                    GetMappabilityRef("hg19")
            
            } else if(input$newref_genome_type == "mm10")  {
                settings_newref$newref_NPA = system.file(
                    "extra-input-files/Mouse_mm10_nonPolyA_ROI.bed", 
                    package = "NxtIRF")
                settings_newref$newref_mappa = 
                    GetMappabilityRef("mm10")
            
            } else if(input$newref_genome_type == "mm9")  {
                settings_newref$newref_NPA = system.file(
                    "extra-input-files/Mouse_mm9_nonPolyA_ROI.bed", 
                    package = "NxtIRF")
                settings_newref$newref_mappa = 
                    GetMappabilityRef("mm9")
            
            } else if(input$newref_genome_type == "(custom)") {
        # do nothing. This allows user to first select the default 
        #   and then change to user-defined files
            } else {
                settings_newref$newref_NPA = ""
                settings_newref$newref_mappa = ""
            }
            settings_newref$ui_newref_genome_type = input$newref_genome_type
        })        
        observeEvent({
            list(
                refresh_tab(), input$release, input$species
            )
        }, {
            req(refresh_tab())
            if(is_valid(input$release) & is_valid(input$species)) {
                withProgress(message = 'Getting options', value = 0, {
                    updateSelectInput(session = session, 
                        inputId = "fasta", 
                        choices = c("", .refresh_genome(
                            input$release, input$species
                        ))
                    )
                    updateSelectInput(session = session, 
                        inputId = "gtf", 
                        choices = c("", .refresh_gtf(
                            input$release, input$species
                        ))
                    )
                })
            } else if(is_valid(input$release)) {
                withProgress(message = 'Getting options', value = 0, {
                    updateSelectInput(session = session, 
                        inputId = "species", 
                        choices = c("", .refresh_species(
                            input$release
                        ))
                    )
                })            
            } else {
                withProgress(message = 'Getting options', value = 0, {
                    updateSelectInput(session = session, 
                        inputId = "release", 
                        choices = c("", .refresh_releases())
                    )
                })
            }
        })
        observeEvent(input$fasta, {
            req(input$fasta)
            settings_newref$newref_fasta = paste0(
                # "https://ftp.ensembl.org/pub/",
                "ftp://ftp.ensembl.org/pub/",
                "release-", as.character(isolate(input$release)),
                "/fasta/", isolate(input$species), "/dna/",
                input$fasta
            )
            output$txt_genome <- renderText(settings_newref$newref_fasta)
        })
        observeEvent(input$gtf, {
            req(input$gtf)
            settings_newref$newref_gtf = paste0(
                # "https://ftp.ensembl.org/pub/",
                "ftp://ftp.ensembl.org/pub/",
                "release-", as.character(isolate(input$release)),
                "/gtf/", isolate(input$species), "/",
                input$gtf
            )
            output$txt_gtf <- renderText(settings_newref$newref_gtf)
        })
        observeEvent(input$buildRef, {
            args <- list(reference_path = settings_newref$newref_path, 
                fasta_file = settings_newref$newref_fasta, 
                gtf_file = settings_newref$newref_gtf,
                genome_type = input$newref_genome_type, 
                nonPolyARef = settings_newref$newref_NPA, 
                MappabilityRef = settings_newref$newref_mappa,
                BlacklistRef = settings_newref$newref_bl
            )

            args <- Filter(is_valid, args)
            if(!("reference_path" %in% names(args))) {
                output$refStatus = renderText({ "Reference path not set" })
            } else if(!any(c("fasta_file") %in% names(args))) {
                output$refStatus = renderText({ "Genome not provided" })        
            } else if(!any(c("gtf_file") %in% names(args))) {
                output$refStatus = renderText("Gene annotations not provided")
            } else {        
                args.df = as.data.frame(t(as.data.frame(args)))
                colnames(args.df) = "value"

                withProgress(message = 'Building Reference', value = 0, {
                    do.call(BuildReference, args)
                })
                # If successfully created, load this reference automatically
                if(file.exists(file.path(settings_newref$newref_path, 
                        "settings.Rds"))) {
                    sendSweetAlert(
                        session = session,
                        title = "Reference Build complete!",
                        type = "success"
                    )           
                } else {
                    sendSweetAlert(
                        session = session,
                        title = paste("Reference Build failed.",
                            "An error must have occurred"),
                        type = "error"
                    )               
                }
            }
        })
            
        # clearNewRef Button
        observeEvent(input$clearNewRef, {
            settings_newref <- setreactive_newref()
            output$txt_reference_path <- 
                renderText("Please select reference path")
            output$txt_genome <- renderText("")
            output$txt_gtf <- renderText("")
            output$txt_mappa <- renderText("")
            output$txt_NPA <- renderText("")
            output$txt_bl <- renderText("")
            updateSelectInput(session = session, 
                inputId = "fasta", 
                choices = c("")
            )
            updateSelectInput(session = session, 
                inputId = "gtf", 
                choices = c("")
            )
            updateSelectInput(session = session, 
                inputId = "species", 
                choices = c("")
            )
            updateSelectInput(session = session, 
                inputId = "release", 
                choices = c("")
            )
        })
        
        return(settings_newref)
    })
}

.refresh_releases <- function() {
    test = XML::getHTMLLinks("http://ftp.ensembl.org/pub")
    test = test[grepl("release-", test)]
    test = test[grepl("/", test)]
    int_release = tstrsplit(test, split="-")[[2]]
    int_release = as.integer(sub("/","",int_release))
    int_release = sort(int_release, decreasing = TRUE)
    return(int_release[int_release > 46])
}

.refresh_species <- function(release) {
    test_genome = XML::getHTMLLinks(paste0(
        "http://ftp.ensembl.org/pub/",
        "release-",
        as.character(release),
        "/fasta/"
    ))
    test_gtf = XML::getHTMLLinks(paste0(
        "http://ftp.ensembl.org/pub/",
        "release-",
        as.character(release),
        "/gtf/"
    ))
    species = union(test_genome, test_gtf)
    species = species[!grepl("..", species, fixed = TRUE)]
    species = sub("/","",species)
    if(all(c("homo_sapiens", "mus_musculus") %in% species)) {
        species = unique(c(
            "homo_sapiens",
            "mus_musculus",
            sort(species)
            )
        )        
    } else {
        species = sort(species)
    }
    return(species)
}

.refresh_genome <- function(release, species) {
    test_genome = XML::getHTMLLinks(paste0(
        "http://ftp.ensembl.org/pub/",
        "release-",
        as.character(release),
        "/fasta/",
        species, "/dna/"           
    ))
    test_genome = test_genome[
        grepl("toplevel", test_genome) |
        grepl("primary", test_genome)
    ]
    test_genome
}

.refresh_gtf <- function(release, species) {
    test_gtf = XML::getHTMLLinks(paste0(
        "http://ftp.ensembl.org/pub/",
        "release-",
        as.character(release),
        "/gtf/",
        species        
    ))
    test_gtf = test_gtf[
        grepl("gtf.gz", test_gtf) &
        !grepl("abinitio", test_gtf) &
        !grepl("scaff", test_gtf) &
        !grepl(".chr.", test_gtf)
    ]
    test_gtf
}