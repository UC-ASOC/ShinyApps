#### #### #### #### #### #### #### #### #### ####
# App: GeoMx-Explorer
# Author: Heewon Seo
#### #### #### #### #### #### #### #### #### ####
# Load libraries
library(shiny)
library(yaml)
library(stringr)
library(GeomxTools)
library(ComplexHeatmap)
library(circlize)

#### #### #### #### #### #### #### #### #### ####
# UI
ui <- fluidPage(
        title = "GeoMx-Explorer ShinyApp",
        theme = "./css/style.css",
        tags$head(tags$link(rel = "shortcut icon", href = "imgs/favicon.ico")),
        h2(HTML(" &#8886; GeoMx-Explorer &#8887;")),
        h4(HTML(" # Check the level of expression with your GeoMx object. ")),
        h4(HTML(" &nbsp; ")),
        fluidRow(
                column(
                        width = 5,
                        fileInput(
                                inputId = "geomxsetobj",
                                label = "Choose a file",
                                multiple = FALSE,
                                accept = c(".rds", ".RDS")
                        ),
                        helpText(h5("Upload a geomx object, e.g., *.GeoMx.RDS")),
                        actionButton(
                                inputId = "exampleData", 
                                label = "Example data", 
                                class = "btn-block",
                                style='padding:4px; font-size:80%; width:100px'
                        ),
                        helpText(h5("Load a subset of the Kidney data from the Nanostring Spatial Organ Atlas."))
                ),
                column(
                        width = 5,
                        textInput(
                                inputId = "genes",
                                label = "Gene symbol(s)",
                                value = "IL32, NPHS1, NPHS2, PODXL, SOD2, LDHA, CYP1B1, HPGDS, TPI1, PPIA, CUBN, SLC5A2, HAVCR1"
                        ),
                        helpText(h5(HTML("<br>The first gene will be shown in the box plot when multiple genes provided.")))
                )
        ),        
        mainPanel(
                tabsetPanel(
                        type = "tabs",
                        tabPanel("Boxplot", h5("The level of expression per sample."), plotOutput("boxplot")),
                        tabPanel("Heat map", h5("The level of expression across samples."), plotOutput("heatmap"))
                )
        ),
        tags$footer( HTML("<footer><div class='copyright'>&copy;2024 GeoMx-Explorer | <a href='https://bit.ly/UC_ASOC' target='_blank'>Applied Spatial Omics Centre</a> | All Rights Reserved.</div></footer>"))
) # End of UI

#### #### #### #### #### #### #### #### #### ####
# Server
server <- function(input, output, session) {
        # options(shiny.maxRequestSize=30*1024^2) # 30 MB

        userdata <- reactiveVal(NULL)
        
        observeEvent(input$exampleData, {
                gSet <- readRDS(file.path("data", "SOA_Kidney.GeoMx.RDS"))
                userdata(gSet)
        })

        observeEvent(input$geomxsetobj, {
                req(file.exists(input$geomxsetobj$datapath))
                gSet <- tryCatch(
                        readRDS(input$geomxsetobj$datapath),
                        error = function(e) e
                )
                validate(
                        need(!inherits(gSet, "error"), "Error reading rds file.")
                )
                userdata(gSet)
        })

        generateBoxplot <- function() {
                if (!is.null(userdata())) {
                        req(input$genes)

                        targetGenes <- unlist(lapply(str_split(input$genes, ","), str_replace_all, pattern = fixed(" "), replacement = ""))
                        if (length(targetGenes) > 1) {
                                gene <- targetGenes[1]
                        } else {
                                gene <- targetGenes
                        }

                        gSet <- userdata()
                        annot <- pData(gSet)

                        if (any(names(assayData(gSet)) %in% "vsn")) {
                                exprMat <- assayData(gSet)$vsn
                        } else if (any(names(assayData(gSet)) %in% "q_norm")) {
                                exprMat <- assayData(gSet)$q_norm
                                exprMat <- log10(exprMat + 1)
                        } else if (any(names(assayData(gSet)) %in% "exprs")) {
                                exprMat <- assayData(gSet)$exprs
                                exprMat <- log10(exprMat + 1)
                        }
                        colnames(exprMat) <- annot$Library
                        
                        tag <- gSet$Tag[1]
                        if (file.exists(file.path("data", paste0(tag, ".yaml")))) {
                                config <- yaml::yaml.load_file(file.path("data", paste0(tag, ".yaml")))

                                if (tolower(config$ColorKey) == "patient") {
                                        colOrder <- config$Patient$Name
                                        col <- config$Patient$Color
                                        names(col) <- colOrder
                                        boxcols <- col[annot$Patient]
                                } else if (tolower(config$ColorKey) == "sample") {
                                        colOrder <- config$Sample$Name
                                        col <- config$Sample$Color
                                        names(col) <- colOrder
                                        boxcols <- col[annot$Sample]
                                } else if (tolower(config$ColorKey) == "region") {
                                        colOrder <- config$Region$Name
                                        col <- config$Region$Color
                                        names(col) <- colOrder
                                        boxcols <- col[annot$Region]
                                } else if (tolower(config$ColorKey) == "segment") {
                                        colOrder <- config$Segment$Name
                                        col <- config$Segment$Color
                                        names(col) <- colOrder
                                        boxcols <- col[annot$Segment]
                                }
                        } else {
                                boxcols <- rep("grey80", ncol(exprMat))
                        }

                        geneCnt <- length(which(rownames(exprMat) == gene))
                        orderIdx <- c(1:ncol(exprMat))
                        if (geneCnt > 0) {
                                geneIdx <- which(rownames(exprMat) == gene)
                                orderIdx <- order(exprMat[geneIdx, ])
                        }

                        par(oma=c(6,0,0,0))
                        boxplot(exprMat[,orderIdx], col=boxcols[orderIdx], pch=20, ylab="Expression", las=3)
                        if (geneCnt > 0) {
                                points(exprMat[geneIdx,orderIdx], col="red", pch=20, cex=2)
                                legend("topright", legend=gene, pch=20, col="red", bty="n")
                                if (!is.null(col)) {
                                        legend("topleft", legend=names(col), fill=col, bty="n")
                                }
                        } else {
                                mtext("Check your gene symbol(s)", side=3, line=1)
                        }
                }
        }

        output$boxplot <- renderPlot({ generateBoxplot() }, height = 600, width = 1200)

        generateHeatmap <- function() {
                if (!is.null(userdata())) {
                        req(input$genes)

                        targetGenes <- unlist(lapply(str_split(input$genes, ","), str_replace_all, pattern = fixed(" "), replacement = ""))                

                        gSet <- userdata()
                        annot <- pData(gSet)

                        if (any(names(assayData(gSet)) %in% "vsn")) {
                                exprMat <- assayData(gSet)$vsn
                        } else if (any(names(assayData(gSet)) %in% "q_norm")) {
                                exprMat <- assayData(gSet)$q_norm
                                exprMat <- log10(exprMat + 1)
                        } else if (any(names(assayData(gSet)) %in% "exprs")) {
				exprMat <- assayData(gSet)$exprs
                                exprMat <- log10(exprMat + 1)
			}
                        colnames(exprMat) <- annot$Library

                        tag <- gSet$Tag[1]

                        subMat <- exprMat[which(rownames(exprMat) %in% targetGenes),]
                        colnames(subMat) <- annot$Library

                        if (file.exists(file.path("data", paste0(tag, ".yaml")))) {
                                config <- yaml::yaml.load_file(file.path("data", paste0(tag, ".yaml")))

                                if (tag == "SOA_Kidney") {
                                        slideOrder <- config$Slide$Name
                                        slideCols <- config$Slide$Color
                                        names(slideCols) <- slideOrder
                                        
                                        sampleOrder <- config$Sample$Name
                                        sampleCols <- config$Sample$Color
                                        names(sampleCols) <- sampleOrder
                                        
                                        regionOrder <- config$Region$Name
                                        regionCols <- config$Region$Color
                                        names(regionCols) <- regionOrder
                                        
                                        segmentOrder <- config$Segment$Name
                                        segmentCols <- config$Segment$Color
                                        names(segmentCols) <- segmentOrder

                                        topAnnotation <- ComplexHeatmap::HeatmapAnnotation(
                                                df = data.frame(
                                                        Slide = factor(annot$Slide, levels = slideOrder),
                                                        Sample = factor(annot$Patient, levels = sampleOrder),
                                                        Region = factor(annot$Region, levels = regionOrder),
                                                        Segment = factor(annot$Segment, levels = segmentOrder),
                                                        Area = annot$area,
                                                        Nuclei = annot$nuclei
                                                ),
                                                col = list(
                                                        Sample = sampleCols,
                                                        Region = regionCols,
                                                        Area = colorRamp2(c(min(annot$area), max(annot$area)), c("white", "grey50")),
                                                        Nuclei = colorRamp2(c(min(annot$nuclei), max(annot$nuclei)), c("white", "blue"))
                                                ),
                                                annotation_legend_param = list(
                                                        Area = list(title = "Area (micrometer^2)")
                                                )
                                        )
                                } else if (tag == "UC_Kidney2021") {
                                        sampleOrder <- config$Sample$Name
                                        sampleCols <- config$Sample$Color
                                        names(sampleCols) <- sampleOrder

                                        regionOrder <- config$Region$Name
                                        regionCols <- config$Region$Color
                                        names(regionCols) <- regionOrder

                                        topAnnotation <- ComplexHeatmap::HeatmapAnnotation(
                                                df = data.frame(
                                                        Sample = factor(annot$Sample, levels = sampleOrder),
                                                        Region = factor(annot$Region, levels = regionOrder),
                                                        Area = annot$area
                                                ),
                                                col = list(
                                                        Sample = sampleCols,
                                                        Region = regionCols,
                                                        Area = colorRamp2(c(min(annot$area), max(annot$area)), c("white", "blue"))
                                                ),
                                                annotation_legend_param = list(
                                                        Area = list(title = "Area (micrometer^2)")
                                                )
                                        )
                                } else if (tag == "UC_mBackbone2023") {
                                        modelOrder <- config$Model$Name
                                        modelCols <- config$Model$Color
                                        names(modelCols) <- modelOrder
                                        
                                        ageOrder <- config$Age$Name
                                        ageCols <- config$Age$Color
                                        names(ageCols) <- ageOrder

                                        regionOrder <- config$Region$Name
                                        regionCols <- config$Region$Color
                                        names(regionCols) <- regionOrder

                                        topAnnotation <- ComplexHeatmap::HeatmapAnnotation(
                                                df = data.frame(
                                                        Model = factor(annot$Model, levels = modelOrder),
                                                        Age = factor(annot$Age, levels = ageOrder),
                                                        Region = factor(annot$Region, levels = regionOrder),
                                                        Area = annot$area,
                                                        Nuclei = annot$nuclei
                                                ),
                                                col = list(
                                                        Model = modelCols,
                                                        Age = ageCols,
                                                        Region = regionCols,
                                                        Area = colorRamp2(c(min(annot$area), max(annot$area)), c("white", "grey50")),
                                                        Nuclei = colorRamp2(c(min(annot$nuclei), max(annot$nuclei)), c("white", "blue"))
                                                ),
                                                annotation_legend_param = list(
                                                        Area = list(title = "Area (micrometer^2)")
                                                )
                                        )
                                } else if (tag == "UC_Liver2023") {
                                        slideOrder <- config$Slide$Name
                                        slideCols <- config$Slide$Color
                                        names(slideCols) <- slideOrder
                                        
                                        patientOrder <- config$Patient$Name
                                        patientCols <- config$Patient$Color
                                        names(patientCols) <- patientOrder

                                        sexOrder <- config$Sex$Name
                                        sexCols <- config$Sex$Color
                                        names(sexCols) <- sexOrder

                                        segmentOrder <- config$Segment$Name
                                        segmentCols <- config$Segment$Color
                                        names(segmentCols) <- segmentOrder

                                        topAnnotation <- ComplexHeatmap::HeatmapAnnotation(
                                                df = data.frame(
                                                        Slide = factor(annot$Slide, levels = slideOrder),
                                                        Patient = factor(annot$Patient, levels = patientOrder),
                                                        Sex = factor(annot$Sex, levels = sexOrder),
                                                        Segment = factor(annot$Segment, levels = segmentOrder),
                                                        Area = annot$area,
                                                        Nuclei = annot$nuclei
                                                ),
                                                col = list(
                                                        Slide = slideCols,
                                                        Patient = patientCols,
                                                        Sex = sexCols,
                                                        Segment = segmentCols,
                                                        Area = colorRamp2(c(min(annot$area), max(annot$area)), c("white", "grey50")),
                                                        Nuclei = colorRamp2(c(min(annot$nuclei), max(annot$nuclei)), c("white", "blue"))
                                                ),
                                                annotation_legend_param = list(
                                                        Area = list(title = "Area (micrometer^2)")
                                                )
                                        )
                                }
                        } else {
                                topAnnotation <- ComplexHeatmap::HeatmapAnnotation(
                                        df = data.frame(
                                                Area = annot$area,
                                                Nuclei = annot$nuclei
                                        ),
                                        col = list(
                                                Area = colorRamp2(c(min(annot$area), max(annot$area)), c("white", "grey50")),
                                                Nuclei = colorRamp2(c(min(annot$nuclei), max(annot$nuclei)), c("white", "blue"))
                                        ),
                                        annotation_legend_param = list(
                                                Area = list(title = "Area (micrometer^2)")
                                        )
                                )
                        }

                        set.seed(10)
                        HM <- ComplexHeatmap::Heatmap(
                                subMat,
                                name = "Expression",
                                cluster_rows = TRUE,
                                show_column_names = FALSE,
                                column_title_side = "bottom",
                                cluster_columns = TRUE,
                                column_title = "",
                                row_title_side = "left",
                                row_title = "",
                                row_names_max_width = unit(6, "cm"),
                                show_heatmap_legend = TRUE,
                                top_annotation = topAnnotation,
                                col = circlize::colorRamp2(c(min(subMat), max(subMat)), c("white", "red"))
                        )

                        draw(HM)
                }
        }

        output$heatmap <- renderPlot({ generateHeatmap() }, height = 600, width = 1200)
} # End of Server

#### #### #### #### #### #### #### #### #### ####
# Allez-y!
shinyApp(ui = ui, server = server)
# (@'3(o.o);
