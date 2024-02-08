import {
    allApplications,
    allowedMonitors,
    allowedResults,
    allTemplates,
// @ts-ignore
} from "@exabyte-io/application-flavors.js";

import { Application, ApplicationMixin } from "./application";
import * as context from "./context";
import { ExecutableContextProvider } from "./context/providers";
import { ContextProviderRegistry } from "./context/registry";
import { Executable, ExecutableMixin } from "./executable";
import { Flavor, FlavorMixin } from "./flavor";
import { Template, TemplateMixin } from "./template";
import { getAllApplications, getApplication } from "./tree";

export {
    Application,
    ApplicationMixin,
    Executable,
    ExecutableMixin,
    Flavor,
    FlavorMixin,
    Template,
    TemplateMixin,
    getAllApplications,
    getApplication,
    allApplications,
    allTemplates,
    allowedResults,
    allowedMonitors,
    ExecutableContextProvider,
    ContextProviderRegistry,
    context,
};
