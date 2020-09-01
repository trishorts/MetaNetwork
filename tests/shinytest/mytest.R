app <- ShinyDriver$new("../../", seed = 0)
app$snapshotInit("mytest")

app$uploadFile(dataFile = "Unsupervised_data.csv") # <-- This should be the path to the file, relative to the app's tests/shinytest directory
app$uploadFile(groupsFile = "Unsupervised_groups.csv") # <-- This should be the path to the file, relative to the app's tests/shinytest directory
app$uploadFile(databaseFile = "filtered_Hsapiens.tab") # <-- This should be the path to the file, relative to the app's tests/shinytest directory
app$setInputs(action = "click")
app$snapshot()
