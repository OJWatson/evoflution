
##' Convert List to XML
##'
##' Can convert list or other object to an xml object using xmlNode
##' @title List to XML
##' @param list - list or list of lists pertaining to the xml
##' @param tag - the opening tag required for the xml creation, e.g. "beast"
##' @return xmlNode
##' @export
#--------------------------------------------------------------------------------------------------#
List_To_Xml <- function(item, tag) {


  # just a textnode, or empty node with attributes
  if(typeof(item) != 'list') {
    if (length(item) > 1) {
      xml <- XML::xmlNode(tag)
      for (name in names(item)) {
        XML::xmlAttrs(xml)[[name]] <- item[[name]]
      }
      return(xml)
    } else {
      return(XML::xmlNode(tag, item))
    }
  }

  # create the node
  if (identical(names(item), c("text", ".attrs"))) {
    # special case a node with text and attributes
    xml <- XML::xmlNode(tag, item[['text']])
  } else {
    # node with child nodes
    xml <- XML::xmlNode(tag)
    for(i in 1:length(item)) {
      if (names(item)[i] != ".attrs") {
        xml <- XML::append.xmlNode(xml, List_To_Xml(item[[i]], names(item)[i]))
      }
    }
  }

  # add attributes to node
  attrs <- item[['.attrs']]
  for (name in names(attrs)) {
    XML::xmlAttrs(xml)[[name]] <- attrs[[name]]
  }
  return(xml)
}
#======================
