## playing with classes 


x <- 1
attr(x, "class") <- "foo"
x

# Or in one line
x <- structure(1, class = "foo")
x


class(x) <- "foo"
class(x)
# [1] "foo"

class(x) <- c("A", "B")
class(x) <- LETTERS

mean <- function (x, ...) {
  UseMethod("mean", x)
}


x <- structure(1, class = letters)
bar <- function(x) UseMethod("bar", x)
bar.z <- function(x) "z"
bar(x)

methods("t")


x <- structure(as.list(1:10), class = "myclass")
length(x)
# [1] 10

mylength.list <- function(x) UseMethod("mylength", x)
mylength <- function(x) length(x)
mylength(x)
# Error in UseMethod("mylength", x) : 
#  no applicable method for 'mylength' applied to an object of class
#  "myclass"


baz <- function(x) UseMethod("baz", x)
baz.D <- function(x) "Z"
baz.B <- function(x) "B"

# so the dot is important for linking up the method to a particular class
baz.F <- function(x) "M"
baz.default <- function(x) "Yo"

ab <- structure(1, class = c("F", "B"))
ba <- structure(1, class = c("O", "O"))
baz(ab)
baz(ba)

# use the method for the next class down
baz.C <- function(x) c("C", NextMethod())
ca <- structure(1, class = c("C", "A"))
cb <- structure(1, class = c("C", "B"))
baz(ca)
baz(cb)


# methods are cool cause they associate with objects of a particular class.
# we can use the "." to sort it out


# will use an S4 class system for repeats
# that way we can check validity and main tain consistency
# it should make the code safer and cleaner

sides <- function(object) 0
setGeneric("sides")


setGeneric("sides", valueClass = "numeric", function(object) {
  standardGeneric("sides")
})

setClass("Shape")
setClass("Polygon", representation(sides = "integer"), contains = "Shape")
setClass("Triangle", contains = "Polygon")
setClass("Square", contains = "Polygon")
setClass("Circle", contains = "Shape")


setMethod("sides", signature(object = "Polygon"), function(object) {
  object@sides
})
setMethod("sides", signature("Triangle"), function(object) 3)
setMethod("sides", signature("Square"),   function(object) 4)
setMethod("sides", signature("Circle"),   function(object) Inf)





showMethods("sides")
showMethods(class = "Triangle")



##### Object oriented programing


setClass("Vehicle")
setClass("Truck", contains = "Vehicle")
setClass("Car", contains = "Vehicle")

setClass("Inspector", representation(name = "character"))
setClass("StateInspector", contains = "Inspector")


setGeneric("inspect.vehicle", function(v, i) {
  standardGeneric("inspect.vehicle")
})


setMethod("inspect.vehicle", 
          signature(v = "Vehicle", i = "Inspector"),  # signiture is the class the method is passed to 
          function(v, i) {
            message("Looking for rust")
          })

setMethod("inspect.vehicle", 
          signature(v = "Car", i = "Inspector"),
          function(v, i) {  
            callNextMethod() # perform vehicle inspection
            message("Checking seat belts")
          })

inspect.vehicle(new("Car"), new("Inspector"))




setMethod("inspect.vehicle", 
          signature(v = "Truck", i = "Inspector"),
          function(v, i) {
            callNextMethod() # perform vehicle inspection
            message("Checking cargo attachments")
          })

inspect.vehicle(new("Truck"), new("Inspector"))
# Looking for rust
# Checking cargo attachments

setMethod("inspect.vehicle", 
          signature(v = "Car", i = "StateInspector"),
          function(v, i) {
            callNextMethod() # perform car inspection
            message("Checking insurance")
          })

inspect.vehicle(new("Car"), new("StateInspector"))


###### Disbatch 2

setClass("C", contains = "character")
setClass("B", contains = "C")
setClass("A", contains = "B")

a <- new("A", "a")
b <- new("B", "b")
c <- new("C", "c")

setGeneric("f", function(x, y) standardGeneric("f"))

setMethod("f", signature("C", "C"), function(x, y) "c-c")
setMethod("f", signature("A", "A"), function(x, y) "a-a")

f(c, c)
f(a, a)

f(b, b)
f(a, c)


setClass("BC", contains = c("B", "C"))
bc <- new("BC", "bc")


setMethod("f", signature("B", "C"), function(x, y) "b-c")
setMethod("f", signature("C", "B"), function(x, y) "c-b")
f(b, b)



repGR <- GRanges(seqnames = Rle(rep$genoChr), 
                 ranges = IRanges(start = rep$genoStart, end = rep$genoEnd),
                 strand = Rle(rep$strand),
                )



setClass(Class = "repeats",contains = "rmsk")
setClass(Class = "rmsk",contains = c("genoRange", "repRange"))

setClass(Class = "genome")
setClass(Class = "hg19", contains = "genome")


setClass(Class = "genoRange", representation(gR = "GRanges"))
setClass(Class = "repRange", representation(rR = "IRanges"))

setClass(Class = "rmsk", representation(geno = "genoRange", rep = "repRange"))

genoGR <- new(Class = "genoRange", repGR)
repGR <- new(Class = "repRange", IRanges(start = rep$repStart, end = rep$repStart + 100))
### this might all be possible within the GRanges package

repObject <- new(Class = "rmsk", rep = repGR , geno = genoGR)



### we can use seqLengths to create the binObjects

### we can use the genome class to get the bin information




