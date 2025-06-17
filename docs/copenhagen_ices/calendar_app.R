
ui <- fluidPage(
  titlePanel("Conference Schedule"),
  fluidRow(
    column(width = 12,
           selectInput("day_select", "Choose Day:", 
                       choices = c("Day 1", "Day 2", "Day 3"), 
                       selected = "Day 1")
    )
  ),
  fluidRow(
    column(width = 12,
           DTOutput("schedule_table")
    )
  )
)

server <- function(input, output, session) {
  output$schedule_table <- renderDT({
    base_path <- dirname(session$clientData$url_pathname)
    base_path <- ifelse(base_path == "/", "", paste0(base_path, "/"))
    
    day_data <- switch(input$day_select,
                       "Day 1" = data.frame(
                         Time = c("11:00 - 11:30", "11:30 - 12:30", "12:30 - 1:00", 
                                  "1:00 - 2:00", "2:00 - 2:30", "2:20 - 3:00", 
                                  "3:00 - 4:00", "4:00 - 4:15", "4:15 - 5:00"),
                         Topic = c(
                           "Session 1: Topic 1",
                           "Session 2a: Topic 2",
                           as.character(a(href = paste0(base_path, "day1_practical.html"), 
                                          "Practical Session 1")),
                           "Lunch break",
                           "Session 3",
                           "Session 3b",
                           as.character(a(href = paste0(base_path, "day1_practical2.html"),
                                          "Practical Session 2")),
                           "Coffee break",
                           as.character(a(href = paste0(base_path, "day1_practical3.html"),
                                          "Practical Session 3")
                           ),
                           Instructor = c("Presenter X", "Presenter Y", "Workshop Lead", 
                                          "", "Presenter Z", "Presenter W", 
                                          "Workshop Lead", "", ""),
                           stringsAsFactors = FALSE
                         ),
                         "Day 2" = data.frame(
                           Time = c("9:00 - 10:00", "10:00 - 11:00", "11:00 - 11:30"),
                           Topic = c(
                             "Day 2 Keynote",
                             "Advanced Topics",
                             as.character(a(href = paste0(base_path, "day2_workshop.html"),
                                            "Workshop 1"))
                           ),
                           Instructor = c("Presenter A", "Presenter B", "Workshop Lead"),
                           stringsAsFactors = FALSE
                         ),
                         "Day 3" = data.frame(
                           Time = c("10:00 - 12:00", "12:00 - 1:00", "1:00 - 3:00"),
                           Topic = c("Final Presentations", "Award Ceremony", "Closing Reception"),
                           Instructor = c("All Teams", "Committee", ""),
                           stringsAsFactors = FALSE
                         )
                       ))
                       
                       datatable(day_data,
                                 escape = FALSE,
                                 rownames = FALSE,
                                 options = list(
                                   dom = 't',
                                   pageLength = 25
                                 ))
  })
}

shinyApp(ui, server)