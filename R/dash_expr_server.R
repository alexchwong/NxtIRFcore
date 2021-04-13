server_expr <- function(id, refresh_tab, volumes, get_threads_reactive, 
        limited = FALSE) {
    moduleServer(id, function(input, output, session) {
    ns = NS(id)
    settings_expr <- setreactive_expr()
    
    observeEvent(list(refresh_tab(), settings_expr$se), {
        output$se_expr_infobox <- renderUI({
            .infobox_update_se(settings_expr$se, settings_expr$collate_path)
        })
    })
    observe({
        shinyDirChoose(input, "dir_reference_path_load", 
            roots = volumes(), session = session)
        shinyDirChoose(input, "dir_bam_path_load", 
            roots = volumes(), session = session)
        shinyDirChoose(input, "dir_irf_path_load", 
            roots = volumes(), session = session)
        shinyDirChoose(input, "dir_collate_path_load", 
            roots = volumes(), session = session)
        shinyFileChoose(input, "file_expr_anno_load", 
            roots = volumes(), session = session)
    })
    observeEvent(input$dir_reference_path_load, {
        req(input$dir_reference_path_load)
        settings_expr$ref_path = parseDirPath(volumes(), 
            input$dir_reference_path_load)
    })
    observeEvent(input$dir_bam_path_load, {
        req(input$dir_bam_path_load)
        settings_expr$bam_path = parseDirPath(volumes(), 
            input$dir_bam_path_load)
    })
    observeEvent(input$dir_irf_path_load, {
        req(input$dir_irf_path_load)
        settings_expr$irf_path = parseDirPath(volumes(), 
            input$dir_irf_path_load)
    })
    observeEvent(input$file_expr_anno_load, {
        req(input$file_expr_anno_load)
        file_selected<-parseFilePaths(volumes(), input$file_expr_anno_load)
        settings_expr$anno_file = as.character(file_selected$datapath)
    })
    observeEvent(input$dir_collate_path_load, {
        req(input$dir_collate_path_load)
        settings_expr$collate_path = parseDirPath(volumes(), 
            input$dir_collate_path_load)
    })
    observeEvent(input$clearLoadRef,{
        settings_expr$ref_path = ""
        output <- .server_expr_clear_ref(output)   
    })
    observeEvent(settings_expr$ref_path, {
        req(settings_expr$ref_path)
        settings_expr$ref_path <- .server_expr_check_ref_path(
            settings_expr$ref_path)
        output <- .server_expr_parse_ref_path(settings_expr$ref_path, output)
        if(is_valid(settings_expr$ref_path)) {
            settings_expr$ref_settings <- readRDS(file.path(
                settings_expr$ref_path, "settings.Rds"))
        } else {
            settings_expr$ref_settings <- NULL
            sendSweetAlert(session = session,
                title = "Invalid Reference Path",
                text = paste("settings.Rds was not found.",
                    "Please check the path contains a valid NxtIRF reference"
                ), type = "error")           
        }
    })
    observeEvent(settings_expr$df.files, {
        req(settings_expr$df.files && is(settings_expr$df.files, "data.frame"))
        req("sample" %in% colnames(settings_expr$df.files))
        settings_expr$df.anno <- .server_expr_sync_df(
            settings_expr$df.files, settings_expr$df.anno)
    })
    observeEvent(settings_expr$df.anno, {
        req(settings_expr$df.anno && is(settings_expr$df.anno, "data.frame"))
        req("sample" %in% colnames(settings_expr$df.anno))
        req(settings_expr$df.files)
        settings_expr$df.files <- .server_expr_sync_df(
            settings_expr$df.anno, settings_expr$df.files)
    })
    ## Handsontable auto-updates settings_expr$df on user edit
    observeEvent(input$hot_files_expr,{
        req(input$hot_files_expr)
        settings_expr$df.files = hot_to_r(input$hot_files_expr) 
    })
    observeEvent(input$hot_anno_expr,{
        req(input$hot_anno_expr)
        settings_expr$df.anno = hot_to_r(input$hot_anno_expr)
    })
    output$hot_files_expr <- renderRHandsontable({
        .server_expr_gen_HOT(settings_expr$df.files)
    })
    output$hot_anno_expr <- renderRHandsontable({
        .server_expr_gen_HOT(settings_expr$df.anno)
    })
    observeEvent(settings_expr$bam_path,{
        settings_expr$df.files <- Expr_Load_BAMs(
            settings_expr$df.files, settings_expr$bam_path, session)
        output$bam_expr_infobox <- Expr_BAM_update_status(
            settings_expr$df.files, settings_expr$bam_path)
    })
    observeEvent(input$run_irf_expr,{
        req(input$run_irf_expr)
        settings_expr$selected_rows <- Expr_IRF_initiate_run(input, session, 
            get_threads_reactive(), 
            isolate(reactiveValuesToList(settings_expr)))
    })
    observeEvent(input$irf_confirm, {
        if(input$irf_confirm == FALSE) {
            settings_expr$selected_rows = c()
            return()
        } else {
            Expr_IRF_actually_run(input, session, get_threads_reactive(), 
                isolate(reactiveValuesToList(settings_expr)))
        }
        settings_expr$selected_rows = c()
        settings_expr$df.files = Expr_Load_IRFs(
            settings_expr$df.files, settings_expr$irf_path)
        output <- .server_expr_check_irf_path(settings_expr$df.files, 
            settings_expr$irf_path, output)
    })
    observeEvent(settings_expr$irf_path,{
        settings_expr$df.files = Expr_Load_IRFs(
            settings_expr$df.files, settings_expr$irf_path)
        output <- .server_expr_check_irf_path(settings_expr$df.files, 
            settings_expr$irf_path, output)
    })
    observeEvent(settings_expr$anno_file,{
        req(settings_expr$anno_file)
        settings_expr$df.anno <- Expr_Load_Anno(settings_expr$df.anno,
            settings_expr$df.files, settings_expr$anno_files)
    })    
    observeEvent(settings_expr$collate_path, {
        if(is_valid(settings_expr$collate_path) &&
            file.exists(file.path(settings_expr$collate_path, "colData.Rds"))) {
            colData.Rds = readRDS(
                file.path(settings_expr$collate_path, "colData.Rds"))
            if(all(c("df.anno", "df.files") %in% names(colData.Rds))) {
                settings_expr$df.files = colData.Rds$df.files
                settings_expr$df.anno = colData.Rds$df.anno
                if("bam_path" %in% names(colData.Rds)) {
                    settings_expr$bam_path = colData.Rds$bam_path
                }
                if("irf_path" %in% names(colData.Rds)) {
                    settings_expr$irf_path = colData.Rds$irf_path
                }
            }
        }
        output <- .server_expr_parse_collate_path(
            reactiveValuesToList(settings_expr), output)
    })
    observeEvent(input$save_expr,{
        req(input$save_expr)
        .server_expr_save_expr(reactiveValuesToList(settings_expr), session)
    })
    observeEvent(input$load_expr,{
        req(input$load_expr)
        if(is_valid(settings_expr$collate_path) &&
            file.exists(file.path(settings_expr$collate_path, "colData.Rds"))) {
            colData.Rds = readRDS(file.path(settings_expr$collate_path, 
                "colData.Rds"))
            req_columns = c("df.anno", "df.files")
            if(all(req_columns %in% names(colData.Rds))) {
                settings_expr$df.files = colData.Rds$df.files
                settings_expr$df.anno = colData.Rds$df.anno
                if("bam_path" %in% names(colData.Rds)) {
                    settings_expr$bam_path = colData.Rds$bam_path
                }
                if("irf_path" %in% names(colData.Rds)) {
                    settings_expr$irf_path = colData.Rds$irf_path
                }
                output$se_expr_infobox <- renderUI({
                    ui_infobox_expr(ifelse(is_valid(settings_expr$se), 2, 1),
                        "Ready to Build Experiment")
                })                
            }
        }
    })
    output$newcol_expr <- renderUI({
        textInput(ns("newcolumnname_expr"), "New Column Name", 
            sprintf("newcol%s", 1 + ncol(settings_expr$df.anno))
        )
    })
    observeEvent(input$addcolumn_expr, {
        req(input$addcolumn_expr)
        df <- isolate(settings_expr$df.anno)
        newcolumn <- eval(parse(text=sprintf('%s(nrow(df))', 
            isolate(input$type_newcol_expr))))
        settings_expr$df.anno <- data.table::setnames(
            cbind(df, newcolumn, stringsAsFactors=FALSE), 
            c(names(df), isolate(input$newcolumnname_expr))
        )
    })
    observeEvent(input$removecolumn_expr, {
        req(input$removecolumn_expr)
        DT <- as.data.table(isolate(settings_expr$df.anno))
        if(isolate(input$newcolumnname_expr) %in% colnames(DT)) {
            message("removing column")
            DT[, c(input$newcolumnname_expr) := NULL]
            settings_expr$df.anno = as.data.frame(DT)
        }
    })
    observeEvent(input$run_collate_expr, {
        req(input$run_collate_expr)
        req(settings_expr$df.files)
        Experiment = na.omit(as.data.table(
            settings_expr$df.files[, c("sample", "irf_file", "cov_file")]
        ))
        reference_path = settings_expr$ref_path
        output_path = settings_expr$collate_path
        if(Expr_CollateData_Validate_Vars(session,
                Experiment, reference_path, output_path)) {
            withProgress(message = 'Collating IRFinder output', value = 0, {
                CollateData(Experiment, reference_path, output_path, 
                    n_threads = get_threads_reactive())
            })
            Expr_Update_colData(settings_expr$collate_path, 
                settings_expr$df.anno, settings_expr$df.files, 
                settings_expr$bam_path, settings_expr$irf_path, 
                session, post_CollateData = TRUE)        
        }
    })
    observeEvent(input$clear_expr, {
        settings_expr$bam_path = ""
        settings_expr$irf_path = ""
        settings_expr$anno_file = ""
        settings_expr$collate_path = ""
        settings_expr$df.files = c()
        settings_expr$df.anno = c()
        settings_expr$se = NULL
    })
    observeEvent(input$build_expr, {
        if(is_valid(settings_expr$collate_path) &&
                file.exists(file.path(
                    settings_expr$collate_path, "colData.Rds"))) {
            colData = as.data.table(settings_expr$df.anno)
            settings_expr$se = MakeSE(settings_expr$collate_path, colData)
        }
    })
# End of Server function
    return(settings_expr)
    })
}

.server_expr_clear_ref <- function(output) {
    output$fasta_source_infobox <- renderInfoBox(infoBox(""))
    output$gtf_source_infobox <- renderInfoBox(infoBox(""))
    output$mappa_source_infobox <- renderInfoBox(infoBox(""))
    output$NPA_source_infobox <- renderInfoBox(infoBox(""))
    output$BL_source_infobox <- renderInfoBox(infoBox(""))
    output$txt_reference_path_load <- renderText("")
    output$ref_expr_infobox <- renderUI(ui_infobox_ref(""))
    return(output)
}

.server_expr_check_ref_path <- function(ref_path) {
    if(is_valid(ref_path)) {
        ref_settings_file <- file.path(ref_path, "settings.Rds")
        if(file.exists(ref_settings_file)) {
            ref_settings = readRDS(ref_settings_file)
            if("reference_path" %in% names(ref_settings)) {
                return(ref_path)
            } else {
                return("")
            }
        } else {
            return("")
        }
    } else {
        return("")
    }
}

.server_expr_parse_ref_path <- function(ref_path, output) {
    if(is_valid(ref_path)) {
        ref_settings_file <- file.path(ref_path, "settings.Rds")
        ref_settings = readRDS(ref_settings_file)
        output <- .server_expr_load_ref(ref_settings, 
            output)
        output$txt_reference_path_load <- renderText(
            ref_path)
        output$ref_expr_infobox <- renderUI(ui_infobox_ref(
            ref_settings_file))
    } else {
        output <- .server_expr_clear_ref(output)
    }
    return(output)
}

.server_expr_sync_df <- function(df1, df2) {
    if(!is_valid(df2)) {
        return(data.frame(sample = df1$sample, stringsAsFactors = FALSE))
    } else {
        df2 = df2[df2$sample %in% df1$sample,]
        df2 = df2[match(df2$sample, df1$sample)]
        return(df2)
    }
}

.server_expr_gen_HOT <- function(df) {
    if(is_valid(df) && is(df, "data.frame")) {
        rhandsontable(df, useTypes = TRUE, stretchH = "all")
    } else {
        NULL
    }
}

.server_expr_load_ref = function(ref_settings, output) {
    ah <- ah_genome_record <- ah_gtf_record <- NULL
    fasta <- gtf <- mappa <- nonPA <- Black <- NULL
    if("ah_genome" %in% names(ref_settings) &&
            is_valid(ref_settings[["ah_genome"]])) {
        ah = AnnotationHub()
        ah_genome_record <- tryCatch({
            basename(ah$sourceurl[
                which(names(ah) == ref_settings[["ah_genome"]])])
        }, error = function(e) NULL)
    }
    if("ah_transcriptome" %in% names(ref_settings) &&
            is_valid(ref_settings[["ah_transcriptome"]])) {
        if(is.null(ah)) ah = AnnotationHub()
        ah_gtf_record <- tryCatch({
            basename(ah$sourceurl[
                which(names(ah) == ref_settings[["ah_transcriptome"]])])
        }, error = function(e) NULL)
    }
    if(is.null(ah_genome_record) && "fasta_file" %in% names(ref_settings)) {
        fasta = basename(ref_settings[["fasta_file"]])
    }
    if(is.null(ah_gtf_record) && "gtf_file" %in% names(ref_settings)) {
        gtf = basename(ref_settings[["gtf_file"]])
    }
    if("MappabilityRef" %in% names(ref_settings)) {
        mappa = basename(ref_settings[["MappabilityRef"]])
    }
    if("nonPolyARef" %in% names(ref_settings)) {
        nonPA = basename(ref_settings[["nonPolyARef"]])
    }
    if("BlacklistRef" %in% names(ref_settings)) {
        Black = basename(ref_settings[["BlacklistRef"]])
    }   
    output <- .server_expr_load_ref_genome(output, ah_genome_record, fasta)
    output <- .server_expr_load_ref_gtf(output, ah_gtf_record, gtf)
    output <- .server_expr_load_ref_misc(output, mappa, nonPA, Black)
    return(output)
}

.server_expr_load_ref_genome <- function(output, ah_genome_record, fasta) {
    if(is_valid(ah_genome_record)) {
        output$fasta_source_infobox <- renderInfoBox({
            infoBox("Genome - AnnotationHub", "", ah_genome_record,
                icon = icon("dna", lib = "font-awesome"),
                color = "green")
        })      
    } else if(is_valid(fasta)) {
        output$fasta_source_infobox <- renderInfoBox({
            infoBox("Genome - User FASTA", "",
                fasta, 
                icon = icon("dna", lib = "font-awesome"),
                color = "green")
        })
    } else {
        output$fasta_source_infobox <- renderInfoBox({
            infoBox("Genome - INVALID", "",
                "", 
                icon = icon("dna", lib = "font-awesome"),
                color = "red")
        })
    }
    return(output)
}

.server_expr_load_ref_gtf <- function(output, ah_gtf_record, gtf) {
    if(is_valid(ah_gtf_record)) {
        output$gtf_source_infobox <- renderInfoBox({
            infoBox("Gene Annotation - AnnotationHub",  "", ah_gtf_record,
                icon = icon("book-medical", lib = "font-awesome"),
                color = "orange")
        })               
    } else if(is_valid(gtf)) {
        output$gtf_source_infobox <- renderInfoBox({
            infoBox("Gene Annotation - User GTF",  "",
                gtf, 
                icon = icon("book-medical", lib = "font-awesome"),
                color = "orange")
        })
    } else {
        output$gtf_source_infobox <- renderInfoBox({
            infoBox("Gene Annotation - INVALID", "",
                "", 
                icon = icon("book-medical", lib = "font-awesome"),
                color = "red")
        })
    }
    return(output)
}

.server_expr_load_ref_misc <- function(output, mappa, nonPA, Black) {
    if(is_valid(mappa)) {
        output$mappa_source_infobox <- renderInfoBox({
            infoBox("Mappability", "",
                mappa, 
                icon = icon("map", lib = "font-awesome"),
                color = "blue")
        })
    } else {
        output$mappa_source_infobox <- renderInfoBox({
            infoBox("Mappability", "",
                "NOT USED", 
                icon = icon("map", lib = "font-awesome"),
                color = "blue")
        })  
    }
    if(is_valid(nonPA)) {
        output$NPA_source_infobox <- renderInfoBox({
            infoBox("Non-PolyA", "",
                nonPA, 
                icon = icon("font", lib = "font-awesome"),
                color = "purple")
        })
    } else {
        output$NPA_source_infobox <- renderInfoBox({
            infoBox("Non-PolyA", "",
                "NOT USED", 
                icon = icon("font", lib = "font-awesome"),
                color = "purple")
        })
    }
    if(is_valid(Black)) {
        output$BL_source_infobox <- renderInfoBox({
            infoBox("BlackList", "",
                Black, 
                icon = icon("list-alt", lib = "font-awesome"),
                color = "red")
        })
    } else {
        output$BL_source_infobox <- renderInfoBox({
            infoBox("BlackList", "",
                "NOT USED", 
                icon = icon("list-alt", lib = "font-awesome"),
                color = "red")
        })
    }
    return(output)
}

Expr_Load_BAMs = function(df.files, bam_path, session) {
# First assume bams are named by subdirectory names
    if(!is_valid(bam_path)) return(df.files)
    temp.DT = FindSamples(bam_path, suffix = ".bam", use_subdir = TRUE)
    if(!is.null(temp.DT) && nrow(temp.DT) > 0) {
        temp.DT = as.data.table(temp.DT)
        if(length(unique(temp.DT$sample)) == nrow(temp.DT)) {
            # Assume subdirectory names designate sample names
        } else {
            temp.DT = as.data.table(FindSamples(
                bam_path, suffix = ".bam", use_subdir = FALSE))
            if(length(unique(temp.DT$sample)) == nrow(temp.DT)) {
            # Else assume bam names designate sample names
            } else {
                sendSweetAlert(session = session,
                    title = "Incompatible BAM file names",
                    text = paste("Could not determine sample names.",
                        "Please ensure either BAMs are uniquely named by",
                        "sample name,",
                        "or its parent directories are uniquely named."
                    ), type = "error")
                temp.DT = NULL
            }
        }
    } else {
        sendSweetAlert(session = session, 
            title = "No BAM files found",
            text = "No BAM files found", type = "error")            
        temp.DT = NULL
    }
    # compile experiment df with bam paths
    if(!is.null(temp.DT) && nrow(temp.DT) > 0)  {
        colnames(temp.DT)[2] = "bam_file"
        if(is_valid(df.files)) {
            df.files = update_data_frame(df.files, temp.DT)
        } else {
            DT = data.table(sample = temp.DT$sample,
                bam_file = "", irf_file = "", cov_file = "")
            DT[temp.DT, on = "sample", c("bam_file") := get("i.bam_file")]
            df.files = as.data.frame(DT)
        }
        return(df.files)
    } else {
        return(df.files)
    }
}

Expr_BAM_update_status <- function(df.files, bam_path, collate_path) {
    if(is_valid(df.files)) {
        if(is_valid(bam_path) &&
                "bam_file" %in% colnames(df.files) && 
                all(file.exists(df.files$bam_file))) {
            return(renderUI(ui_infobox_bam(bam_path, 
                    df.files$bam_file)))
        } else if("irf_file" %in% colnames(df.files) && 
                all(file.exists(df.files$irf_file))) {
            return(renderUI(ui_infobox_bam(bam_path, escape = TRUE)))
        } else if(is_valid(collate_path) && 
                file.exists(file.path(collate_path, "colData.Rds"))
        ){
            return(renderUI(ui_infobox_bam(bam_path, escape = TRUE)))
        } else if("bam_file" %in% colnames(df.files)) {
            return(renderUI(ui_infobox_bam(bam_path, df.files$bam_file)))
        }
    } else {
        return(renderUI(ui_infobox_bam(bam_path)))        
    } 
}

Expr_Load_IRFs = function(df.files, irf_path) {
    if(!is_valid(irf_path)) return(df.files)
    # merge irfinder paths
    temp.DT = FindSamples(irf_path, suffix = ".txt.gz", use_subdir = FALSE)
    if(!is.null(temp.DT) && nrow(temp.DT) > 0) {
        temp.DT = as.data.table(temp.DT)
        if(length(unique(temp.DT$sample)) == nrow(temp.DT)) {
            # Assume output names designate sample names
        } else {
            temp.DT = as.data.table(FindSamples(
                irf_path, suffix = ".txt.gz", use_subdir = TRUE))
            if(length(unique(temp.DT$sample)) == nrow(temp.DT)) {
            # Else assume subdirectory names designate sample names
            } else {
                temp.DT = NULL
            }
        }
    } else {
        temp.DT = NULL
    }
    if(!is.null(temp.DT) && nrow(temp.DT) > 0) {
        colnames(temp.DT)[2] = "irf_file"
        if(is_valid(df.files)) {
            df.files = update_data_frame(df.files, temp.DT)
        } else {
            DT = data.table(sample = temp.DT$sample,
                bam_file = "", irf_file = "", cov_file = "")
            DT[temp.DT, on = "sample", c("irf_file") := get("i.irf_file")] 
            df.files = as.data.frame(DT)      
        }   
    }
    temp.DT2 = FindSamples(irf_path, suffix = ".cov", use_subdir = FALSE)
    if(!is.null(temp.DT2) && nrow(temp.DT2) > 0) {
        temp.DT2 = as.data.table(temp.DT2)
        if(length(unique(temp.DT2$sample)) == nrow(temp.DT2)) {
            # Assume output names designate sample names
        } else {
            temp.DT2 = as.data.table(FindSamples(
                irf_path, suffix = ".cov", use_subdir = TRUE))
            if(length(unique(temp.DT2$sample)) == nrow(temp.DT2)) {
        # Else assume subdirectory names designate sample names
            } else {
                temp.DT2 = NULL
            }
        }
    } else {
        temp.DT2 = NULL
    }
# compile experiment df with irfinder paths
    if(!is.null(temp.DT2) && nrow(temp.DT2) > 0) {
        colnames(temp.DT2)[2] = "cov_file"
        df.files = update_data_frame(df.files, temp.DT2)
    }
    return(df.files)
}

Expr_IRF_initiate_run <- function(input, session, n_threads, settings_expr) {
    if(!is_valid(settings_expr$df.files)) {
        sendSweetAlert(session = session, type = "error",
            title = "No bam files in experiment",
            text = "Please select bam folder and select bam files")
        return()
    }
    if(!("bam_file" %in% colnames(settings_expr$df.files))) {
        sendSweetAlert(session = session, type = "error",
            title = "No bam files in experiment",
            text = "Please select bam folder and select bam files")
        return()
    }
    if(!is_valid(input$hot_files_expr_select$select$r)) {
        sendSweetAlert(session = session, type = "error",
            title = "No BAM files selected", 
            text = "Please highlight cells of bam files to run IRFinder")
        return()        
    }
    selected_rows = seq(input$hot_files_expr_select$select$r,
        input$hot_files_expr_select$select$r2)
    selected_cols = seq(input$hot_files_expr_select$select$c,
        input$hot_files_expr_select$select$c2)
    bam_col = which(colnames(settings_expr$df.files) == "bam_file")
    bam_files = settings_expr$df.files$bam_file[selected_rows]
    if(!is_valid(settings_expr$ref_path)) {
        sendSweetAlert(session = session,
            title = "Missing Reference", type = "error",
            text = "Please load Reference before running IRFinder")
    } else if(!(bam_col %in% selected_cols)) {
        sendSweetAlert(session = session, type = "error",
            title = "No BAM files selected",
            text = "Please highlight cells of bam files to run IRFinder")
    } else if(!all(file.exists(bam_files))) {
        sendSweetAlert(session = session,
            title = "Missing BAMs", type = "error",
            text = "Please check all selected bam files exist")          
    } else if(!file.exists(file.path(
            settings_expr$ref_path, "IRFinder.ref.gz"))) {
        sendSweetAlert(session = session, type = "error",
            title = "Missing IRFinder Reference",
            text = "IRFinder.ref.gz is missing")
    } else if(!is_valid(settings_expr$irf_path) || 
            !dir.exists(settings_expr$irf_path)) {
        sendSweetAlert(session = session, type = "error",
            title = "Missing IRFinder output path",
            text = "Please set IRFinder output path")
    } else {
        n_threads = min(n_threads, length(selected_rows))
        if(n_threads < length(selected_rows)) {
            n_rounds = ceiling(length(selected_rows) / n_threads)
            n_threads = ceiling(length(selected_rows) / n_rounds)
        }
        msg = paste("Run IRFinder on", length(selected_rows), "samples?",
            "Estimated runtime", 10 * 
                ceiling(length(selected_rows) / n_threads),
            "minutes using", n_threads, 
            "threads (10min per BAM @ 100 million reads per sample)"
        )
        ask_confirmation(inputId = "irf_confirm", type = "warning", 
            title = msg, btn_labels = c("Cancel", "Run IRFinder"),
            btn_colors = c("#00BFFF", "#FE2E2E"))
        return(selected_rows)
    }
}

Expr_IRF_actually_run <- function(input, session, n_threads, settings_expr) {
    n_threads = min(n_threads, length(settings_expr$selected_rows))
    n_rounds = ceiling(length(settings_expr$selected_rows) / n_threads)
    n_threads = ceiling(length(settings_expr$selected_rows) / n_rounds)
    if(n_threads == 1) {
        # run IRFinder using single thread
        withProgress(message = 'Running IRFinder', value = 0, {
            i_done = 0
            incProgress(0.001, 
                message = paste('Running IRFinder',
                    i_done, "of", length(settings_expr$selected_rows), "done")
            )
            for(i in settings_expr$selected_rows) {
                IRFinder(
                    bamfiles = settings_expr$df.files$bam_file[i],
                    sample_names = settings_expr$df.files$sample[i],
                    reference_path = settings_expr$ref_path,
                    output_path = settings_expr$irf_path,
                    n_threads = 1,
                    run_featureCounts = FALSE            
                )
                i_done = i_done + 1
                incProgress(1 / length(settings_expr$selected_rows), 
                    message = paste(i_done, "of", 
                        length(settings_expr$selected_rows), "done")
                )
            }
        })
    } else if(n_threads <= length(settings_expr$selected_rows)) {
        # extract subset to run in parallel
        row_starts = seq(settings_expr$selected_rows[1], by = n_threads,
            length.out = n_rounds)
        withProgress(message = 'Running IRFinder - Multi-threaded', value = 0, {
            i_done = 0
            incProgress(0.001, 
                message = paste('Running IRFinder - Multi-threaded,',
                i_done, "of", length(settings_expr$selected_rows), "done")
            )
            for(i in seq_len(n_rounds)) {
                selected_rows_subset = seq(row_starts[i], 
                    min(length(settings_expr$selected_rows), 
                    row_starts[i] + n_threads - 1)
                )
                IRFinder(
                    bamfiles = settings_expr$df.files$
                        bam_file[selected_rows_subset],
                    sample_names = settings_expr$df.files$
                        sample[selected_rows_subset],
                    reference_path = settings_expr$ref_path,
                    output_path = settings_expr$irf_path,
                    n_threads = n_threads,
                    run_featureCounts = FALSE                  
                )                        
                i_done = i_done + n_threads
                incProgress(n_threads / length(settings_expr$selected_rows), 
                    message = paste(i_done, "of", 
                        length(settings_expr$selected_rows), "done")
                )
            }
        })
    }
    sendSweetAlert(
        session = session,
        title = "IRFinder run completed",
        type = "success"
    )
}

.server_expr_check_irf_path <- function(df.files, irf_path, output) {
    if(is_valid(df.files) && "irf_file" %in% colnames(df.files)) {
        irf_files = df.files$irf_file
    } else {
        irf_files = NULL
    }
    output$irf_expr_infobox <- renderUI({
        ui_infobox_irf(irf_path, irf_files)
    })
    return(output)
}

Expr_Load_Anno = function(df.anno, df.files, anno_file, session) {
    temp.DT = tryCatch(fread(anno_file), error = function(e) NULL)
    if(!is_valid(temp.DT)) return(df.anno)
    if(nrow(temp.DT) == 0) return(df.anno)
    if(!("sample" %in% colnames(temp.DT))) {
        sendSweetAlert(
            session = session,
            title = "Error in Annotation file",
            text = "'sample' must be the name of the first column",
            type = "error"
        )
        return(df.anno)
    }
    files_header = c("bam_file", "irf_file", "cov_file")
    anno_header = names(temp.DT)[!(names(temp.DT) %in% files_header)]
    temp.DT.files = copy(temp.DT)
    if(length(anno_header) > 0) temp.DT.files[, c(anno_header) := NULL]
    if(is_valid(df.files)) {
        df.files = update_data_frame(df.files, temp.DT.files)
    } else {
        DT = data.table(sample = temp.DT$sample, bam_file = "", irf_file = "",
            cov_file = "")
        df.files = update_data_frame(DT, temp.DT.files)
    }
    temp.DT.anno = copy(temp.DT)
    files_header_exist = intersect(files_header, names(temp.DT))
    if(length(files_header_exist) > 0) {
        temp.DT.anno[, c(files_header_exist):= NULL]
    }
    if(is_valid(df.anno)) {
        df.anno = update_data_frame(df.anno, temp.DT.anno)
    } else {
        df.anno = temp.DT.files
    }
    return(df.anno)
}

.server_expr_parse_collate_path <- function(settings_expr, output) {
    if(is_valid(settings_expr$collate_path) &&
        file.exists(file.path(settings_expr$collate_path, "colData.Rds"))) {
        colData.Rds = readRDS(
            file.path(settings_expr$collate_path, "colData.Rds"))
        output$se_expr_infobox <- renderUI({
            ui_infobox_expr(ifelse(is_valid(settings_expr$se),2,1),
                "Ready to Build Experiment")
        })
    } else if(is_valid(settings_expr$collate_path) &&
            is_valid(settings_expr$df.files) &&
            all(file.exists(settings_expr$df.files$irf_file))) {
        output$se_expr_infobox <- renderUI({
            ui_infobox_expr(1, "Ready to run NxtIRF-Collate")
        })
    } else if(is_valid(settings_expr$collate_path)) {
        output$se_expr_infobox <- renderUI({
            ui_infobox_expr(1, "IRFinder files incomplete")
        })
    } else {
        output$se_expr_infobox <- renderUI(ui_infobox_expr(0))
    }
    return(output)
}

.server_expr_save_expr <- function(settings_expr, session) {
    if(is_valid(settings_expr$collate_path) &&
        file.exists(file.path(settings_expr$collate_path, "colData.Rds"))) {
        colData.Rds = list(
            df.anno = settings_expr$df.anno,
            df.files = settings_expr$df.files,
            bam_path = settings_expr$bam_path,
            irf_path = settings_expr$irf_path
        )
        saveRDS(colData.Rds, file.path(settings_expr$collate_path, 
            "colData.Rds"))
        sendSweetAlert(
            session = session,
            title = paste("Annotations saved to", settings_expr$collate_path),
            type = "success"
        )
    } else {
        sendSweetAlert(
            session = session,
            title = "Annotations not saved; run CollateData first!",
            type = "error"
        )
    }
}

Expr_CollateData_Validate_Vars <- function(session,
        Experiment, reference_path, output_path) {
    if(!is_valid(reference_path)) {
        sendSweetAlert(
            session = session,
            title = "Missing Reference",
            text = "Please load Reference before running NxtIRF::CollateData",
            type = "error"
        )
        return(FALSE)
    } else if(!is_valid(output_path)) {
        sendSweetAlert(
            session = session,
            title = "Missing NxtIRF Path",
            text = paste("Please select NxtIRF path before",
                "running NxtIRF::CollateData"),
            type = "error"
        )
        return(FALSE)
    } else if(!dir.exists(output_path)) {
        sendSweetAlert(
            session = session,
            title = "Invalid NxtIRF Path",
            text = "Please make sure NxtIRF output path exists",
            type = "error"
        )
        return(FALSE)
    } else if(nrow(Experiment) == 0) {
        sendSweetAlert(
            session = session,
            title = "No samples found to collate NxtIRF Experiment",
            text = "Please load IRFinder output of some samples",
            type = "error"
        )
        return(FALSE)
    }
    return(TRUE)
}

Expr_Update_colData <- function(collate_path, df.anno, df.files, 
        bam_path, irf_path, session, post_CollateData = FALSE) {
    if(file.exists(file.path(collate_path, "colData.Rds"))) {
        colData.Rds = readRDS(file.path(collate_path, "colData.Rds"))
        if(all(colData.Rds$df.anno$sample %in% df.anno$sample)) {
            colData.Rds$df.anno = df.anno
            colData.Rds$df.files = df.files
            colData.Rds$bam_path = bam_path
            colData.Rds$irf_path = irf_path
            saveRDS(colData.Rds, file.path(collate_path, "colData.Rds"))
            if(post_CollateData) {
                sendSweetAlert(
                    session = session,
                    title = "NxtIRF-Collate run completed",
                    type = "success"
                )
            }
        } else {
            if(post_CollateData) {
                sendSweetAlert(
                    session = session,
                    title = "NxtIRF-Collate did not collate all samples",
                    type = "warning"
                )
            }
        }
    } 
}

.infobox_update_se <- function(se, path) {
    ui_infobox_expr(ifelse(
        is(se, "NxtSE"), 2, ifelse(
            is_valid(path) && file.exists(file.path(path,"colData.Rds")),
            1,0)))
}