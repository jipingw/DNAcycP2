#' @importFrom basilisk BasiliskEnvironment
env1 <- BasiliskEnvironment("env1", pkgname="DNAcycP2",
                            packages=c(
                                'python==3.11',
                                'pandas==2.1.2',
                                'tensorflow==2.14.0',
                                'keras==2.14.0',
                                'docopt==0.6.2'),
                            pip=c(
                                'numpy==1.26.1',
                                'bio==1.7.1'),
                            path="dnacycp_python")
