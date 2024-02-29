import {
    allowedApplications,
    allowedMonitors,
    allowedPostProcessors,
    allowedResults,
    allTemplates,
} from "@exabyte-io/application-flavors.js";

import { Application } from "./application";
import * as context from "./context";
import { ExecutableContextProvider } from "./context/providers";
import { ContextProviderRegistry } from "./context/registry";
import { Executable } from "./executable";
import { Flavor } from "./flavor";
import { Template } from "./template";
import { getAllApplications, getApplication } from "./tree";

export {
    Application,
    Executable,
    Flavor,
    Template,
    getAllApplications,
    getApplication,
    allTemplates,
    allowedApplications,
    allowedResults,
    allowedMonitors,
    allowedPostProcessors,
    ExecutableContextProvider,
    ContextProviderRegistry,
    context,
};
