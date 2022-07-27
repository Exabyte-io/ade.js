import { allApplications, allTemplates, allowedResults, allowedMonitors } from "@exabyte-io/application-flavors.js";

import { Application } from "./application";
import { Executable } from "./executable";
import { Flavor } from "./flavor";
import { Template } from "./template";
import { getAllApplications, getApplication } from "./tree";
import { ExecutableContextProvider } from "./context/providers";
import * as context from "./context";
import { ContextProviderRegistry } from "./context/registry";

export {
    Application,
    Executable,
    Flavor,
    Template,
    getAllApplications,
    getApplication,
    allApplications,
    allTemplates,
    allowedResults,
    allowedMonitors,
    ExecutableContextProvider,
    ContextProviderRegistry,
    context,
}
